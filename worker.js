// JavaScript port of the optimisation in hedgecalculator.py.
// Runs in a Web Worker so the page stays responsive during long runs.
// Points on the sphere are unit vectors in R^3; P is a flat Float64Array [x0,y0,z0, x1,y1,z1, ...].

const R_EARTH = 6371;

// Geometric helpers
function latlon2vec(lat, lon) {
  const a = lat * Math.PI / 180, b = lon * Math.PI / 180;
  return [Math.cos(a) * Math.cos(b), Math.cos(a) * Math.sin(b), Math.sin(a)];
}
function vec2latlon(v) {
  return [Math.asin(Math.max(-1, Math.min(1, v[2]))) * 180 / Math.PI, Math.atan2(v[1], v[0]) * 180 / Math.PI];
}
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const norm = (a) => Math.sqrt(dot(a, a));
function normalize(a) { const n = norm(a); return [a[0] / n, a[1] / n, a[2] / n]; }
// Project v onto the tangent space at x (also used as vector transport)
function proj(x, v) { const d = dot(x, v); return [v[0] - d * x[0], v[1] - d * x[1], v[2] - d * x[2]]; }
// Exponential map on the sphere
function expmap(x, v) {
  const t = norm(v);
  if (t < 1e-16) return x;
  const c = Math.cos(t), s = Math.sin(t) / t;
  return normalize([c * x[0] + s * v[0], c * x[1] + s * v[1], c * x[2] + s * v[2]]);
}
// Central angle between unit vectors x and p (atan2 form, as in the Python version)
function angle(x, p0, p1, p2) {
  const c0 = p1 * x[2] - p2 * x[1];
  const c1 = p2 * x[0] - p0 * x[2];
  const c2 = p0 * x[1] - p1 * x[0];
  return Math.atan2(Math.sqrt(c0 * c0 + c1 * c1 + c2 * c2), p0 * x[0] + p1 * x[1] + p2 * x[2]);
}

// Cost functions. Each returns {f} or {f, g} with g the Riemannian gradient at x.
// "sum":   mean great-circle distance (km)
// "hedge": -(mean score)/5000 = -1/N * sum exp(-10 d / MAX_DIST)
function makeCost(kind, P, maxDist) {
  const N = P.length / 3;
  const k = 10 * R_EARTH / maxDist;
  return function (x, wantGrad) {
    let f = 0, g0 = 0, g1 = 0, g2 = 0;
    for (let i = 0; i < 3 * N; i += 3) {
      const p0 = P[i], p1 = P[i + 1], p2 = P[i + 2];
      const th = angle(x, p0, p1, p2);
      let w;
      if (kind === "sum") { f += th; w = 1; }
      else { w = Math.exp(-k * th); f += w; w *= k; }
      if (wantGrad) {
        // d(theta)/dx on the sphere = -(p - cos(theta) x) / sin(theta)
        const c = p0 * x[0] + p1 * x[1] + p2 * x[2];
        const u0 = p0 - c * x[0], u1 = p1 - c * x[1], u2 = p2 - c * x[2];
        const n = Math.sqrt(u0 * u0 + u1 * u1 + u2 * u2);
        if (n > 1e-12) { g0 -= w * u0 / n; g1 -= w * u1 / n; g2 -= w * u2 / n; }
      }
    }
    if (kind === "sum") {
      const s = R_EARTH / N;
      return { f: s * f, g: [s * g0, s * g1, s * g2] };
    }
    return { f: -f / N, g: [g0 / N, g1 / N, g2 / N] };
  };
}

// Riemannian conjugate gradient (Polak-Ribiere+, Armijo backtracking)
function conjugateGradient(cost, x0, maxIter = 1000, minGrad = 1e-6) {
  let x = x0;
  let { f, g } = cost(x, true);
  let d = [-g[0], -g[1], -g[2]];
  let gg = dot(g, g);
  let lastStep = null;
  for (let it = 0; it < maxIter; it++) {
    if (Math.sqrt(gg) < minGrad) break;
    let slope = dot(g, d);
    if (slope >= 0) { d = [-g[0], -g[1], -g[2]]; slope = -gg; }
    const dn = norm(d);
    // Step length in radians along d: start from twice the last accepted step, capped at a quarter turn
    let alpha = Math.min(lastStep === null ? 0.1 : 2 * lastStep, Math.PI / 2) / dn;
    let xn = null, fn = Infinity;
    for (let bt = 0; bt < 50; bt++) {
      const cand = expmap(x, [alpha * d[0], alpha * d[1], alpha * d[2]]);
      const fc = cost(cand, false).f;
      if (fc <= f + 1e-4 * alpha * slope) { xn = cand; fn = fc; break; }
      alpha *= 0.5;
    }
    if (xn === null || alpha * dn < 1e-13) break;
    lastStep = alpha * dn;
    const r = cost(xn, true);
    const gOld = proj(xn, g), dOld = proj(xn, d);
    const ggNew = dot(r.g, r.g);
    const beta = Math.max(0, (ggNew - dot(r.g, gOld)) / gg);
    d = [-r.g[0] + beta * dOld[0], -r.g[1] + beta * dOld[1], -r.g[2] + beta * dOld[2]];
    x = xn; f = r.f; g = r.g; gg = ggNew;
  }
  return { point: x, cost: f };
}

// Karcher mean: minimiser of the sum of squared geodesic distances
// (what pymanopt's compute_centroid approximates)
function karcherMean(P) {
  const N = P.length / 3;
  let s = [0, 0, 0];
  for (let i = 0; i < 3 * N; i += 3) { s[0] += P[i]; s[1] += P[i + 1]; s[2] += P[i + 2]; }
  let x = norm(s) > 1e-9 ? normalize(s) : [P[0], P[1], P[2]];
  for (let it = 0; it < 200; it++) {
    let v0 = 0, v1 = 0, v2 = 0;
    for (let i = 0; i < 3 * N; i += 3) {
      const p0 = P[i], p1 = P[i + 1], p2 = P[i + 2];
      const c = p0 * x[0] + p1 * x[1] + p2 * x[2];
      const u0 = p0 - c * x[0], u1 = p1 - c * x[1], u2 = p2 - c * x[2];
      const n = Math.sqrt(u0 * u0 + u1 * u1 + u2 * u2);
      if (n > 1e-12) {
        const th = Math.atan2(n, c) / n; // log map: theta * u / |u|
        v0 += th * u0; v1 += th * u1; v2 += th * u2;
      }
    }
    const v = [v0 / N, v1 / N, v2 / N];
    x = expmap(x, v);
    if (norm(v) < 1e-12) break;
  }
  return x;
}

// Run CG from n_seeds randomly chosen data points and keep the best.
// Seeds are drawn with replacement like np.random.choice; duplicates give identical runs, so they are skipped.
function multiStart(cost, P, nSeeds, onProgress) {
  const N = P.length / 3;
  const seeds = new Set();
  for (let i = 0; i < nSeeds; i++) seeds.add(Math.floor(Math.random() * N));
  let best = null, done = 0;
  for (const s of seeds) {
    const res = conjugateGradient(cost, [P[3 * s], P[3 * s + 1], P[3 * s + 2]]);
    if (best === null || res.cost < best.cost) best = res;
    onProgress(++done, seeds.size);
  }
  return best;
}

self.onmessage = (e) => {
  const { coords, maxDist, nSeeds } = e.data;
  const N = coords.length;
  const P = new Float64Array(3 * N);
  coords.forEach(([lat, lon], i) => P.set(latlon2vec(lat, lon), 3 * i));

  self.postMessage({ type: "progress", text: "Computing centroid…" });
  const centroid = karcherMean(P);

  const sumCost = makeCost("sum", P, maxDist);
  const avgPoint = multiStart(sumCost, P, nSeeds, (i, n) =>
    self.postMessage({ type: "progress", text: `Min sum of distances: seed ${i}/${n}` })).point;

  const hedgeCost = makeCost("hedge", P, maxDist);
  const hedge = multiStart(hedgeCost, P, nSeeds, (i, n) =>
    self.postMessage({ type: "progress", text: `Best hedge: seed ${i}/${n}` })).point;

  const score = (x) => -5000 * hedgeCost(x, false).f;
  const meanDist = (x) => sumCost(x, false).f;
  const describe = (x) => ({ latlon: vec2latlon(x), score: score(x), meanDist: meanDist(x) });
  self.postMessage({
    type: "result",
    hedge: describe(hedge),
    centroid: describe(centroid),
    avgPoint: describe(avgPoint),
  });
};
