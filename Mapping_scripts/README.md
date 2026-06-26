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

--------------------------------------------------
## E3SM to GCAM land mapping
### Script name - `e3sm_gcam_land_mapping.r`

* Generate glm input files and iESMv2 IAC grid mapping files that reflect current GCAM Regions/GLUs
* For iESM in E3SM v2 and v3, using the original GLM
* Both the region and USA-state `elm2gcam` files are written

--------------------------------------------------
## E3SM to GCAM land mapping
### Script name - `elm2gcam_mapping_generator.py`

Python translation of the `elm2gcam` portion of `e3sm_gcam_land_mapping.r`.

--------------------------------------------------
## ERA5 HDD/CDD to GCAM land mapping
### Script name - `era5_hdd_cdd_gcam_regions.py`

Maps ERA5 mean-annual HDD and CDD to GCAM land units (region x GLU),
including a second set of output files that sub-divide the USA into
individual states following the same approach as `e3sm_gcam_land_mapping.r`.

#### APPROACH:
The MOIRAI high-resolution land-area raster (5 arcminute, 2160x4320)
encodes each pixel's GCAM land unit as region_code*10000 + GLU_code.
ERA5 HDD/CDD fields are interpolated to that same fine grid, then a
land-area-weighted mean is computed for every GCAM land unit.
A second raster (moirai_state_rast) replaces the USA region code inside
the contiguous US with individual state codes, enabling state-level output.

--------------------------------------------------
## ERA5 HDD/CDD to GCAM xml file creaton
### Script name - `era5_hdd_cdd_to_gcam_xml.py`

Converts ERA5-derived HDD and CDD CSVs (output of `era5_hdd_cdd_gcam_regions.py`)
into GCAM XML input files that follow the same structure as:
   `HDDCDD_constdd_no_GCM.xml`   (32 global GCAM regions)
   `HDDCDD_constdds_USA.xml`     (51 US states)

APPROACH:
The ERA5 CSVs contain one row per region x GLU (or state x GLU).
These are collapsed to one value per region (or state) using a
land-area-weighted mean across GLUs.  The result is rounded to the
nearest integer to match the existing XML format.

ERA5 PERIOD LABEL → XML HISTORICAL YEAR MAPPING:
  ERA5 1980 (avg 1976–1980)  → XML year 1975
  ERA5 1990 (avg 1986–1990)  → XML year 1990
  ERA5 2005 (avg 2001–2005)  → XML year 2005
  ERA5 2010 (avg 2006–2010)  → XML year 2010
  ERA5 2015 (avg 2011–2015)  → XML year 2015
  Future years 2020–2100     → hold 2015 value constant