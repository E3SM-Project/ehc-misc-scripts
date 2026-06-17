# Scripts for generating mapping files
--------------------------------------------------
## Country to grid mapping
### Script name - `country2grid_mapping.r`

An R utility script designed for spatial downscaling, grid mapping, and region aggregation between the **E3SM Land Model (ELM)** and the **Global Change Analysis Model (GCAM)**. 

The script generates grid mapping files that match current GCAM regions and Geographic Land Units (GLUs) to half-degree and nominal 1-degree resolution grids used by the Global Land Model (GLM) and ELM.


### Features

* **Grid Downscaling & Weighting:** Extracts spatial overlap data, computing the exact area fractions (`Weight`) of fine-resolution country/GLU data falling within coarser climate model grid cells.
* **Antimeridian Handling:** Robustly manages grid cells that straddle the $+180^{\circ}/-180^{\circ}$ longitude boundary by dynamically splitting them into dual spatial objects.
* **Multi-Core Parallelization:** Leverages the `parallel` library (`mclapply`) to accelerate intensive spatial grid intersections across all available CPU cores.
* **Format Diversity:** Processes a mix of spatial data formats, including binary raster grids (`.bil`), ESRI shapefile lookups (`.csv`), and NetCDF parameters (`.nc`).


### Output files generated
| File Name | Description |
|----------|-------------------|
| `elm0.9x1.25tocountry_mapping.csv` | The primary grid mapping of GCAM regions/GLUs to the ELM nominal 1-degree grid with calculated cell weights. | 
| `gcam_region_grid.txt` / `.asc` | Half-degree GCAM region map as raw text (no header) and ESRI ASCII format. | 
| `gcam_zone_grid.txt` | Half-degree GCAM GLU map enumerated within each region. | 
| `gcam2glm_mapping.csv` | Grid mapping of GCAM region/GLU combos down to the GLM half-degree grid with weights. |
 
 
### Usage
```
Rscript country2grid_mapping.r
```
---