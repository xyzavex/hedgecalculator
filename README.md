# Hedge Calculator

A small website for finding the best place to "hedge" in GeoGuessr. If you know the country but have no other clues, where should you click to get the highest average score?

Upload a JSON file of location and the site calculates three points:

- **Best Hedge**: maximizes your average score over the locations.
- **Min Sum of Dists**: minimizes the average distance to the locations.
- **Centroid**: minimizes the average squared distance to the locations.

It shows their coordinates and plots them on a map with the locations. You can also change two settings: `MAX_DIST`, the distance scale used in the scoring formula (default 17502 km, estimated from the map "A Competitive World"), and `n_seeds`, the number of random starting points for the optimizer (default 100).

Everything runs in your browser, so your file is never uploaded anywhere.

## Pages

- **Calculator** (`index.html`): upload a file and run the calculation.
- **Docs** (`docs.html`): description of the site. 
- **Examples** (`examples.html`): precomputed results for 13 countries and 3 water hedges.

## Files

| File | Purpose |
|---|---|
| `index.html` | Calculator page: file upload, settings, results table and map drawing |
| `worker.js` | The optimization, which runs in a background thread so the page stays responsive |
| `docs.html` | Docs page |
| `examples.html` | Examples page |
| `examples/` | Figures used on the Examples page |
| `style.css` | Shared styles and navigation bar |

The map uses [d3](https://d3js.org/) and [Natural Earth](https://www.naturalearthdata.com/) outlines from [world-atlas](https://github.com/topojson/world-atlas). All three are loaded from the jsDelivr CDN, so there is nothing to install.

## About the code

The code for this was originally written in Python, using pymanopt for the optimization and cartopy for the maps. I used AI to convert it into a website and change the implementation to Javascript. Its results were checked against the Python version and match to six decimal places.
