# Scripts for regridding input data
--------------------------------------------------

## Create gridded CO2 emissision files for E3SM-GCAM simulation:
### Script name - `generate_initial_co2_files.sh`

* Create baseline (2014) gridded CO2 emission files for E3SM-GCAM
* One input argument determines resolution
* Three main output csv files: aircraft, shipment, surface
* Also generates associated netcdf files and some intermediate files


## Create population density file at various resolutions:
### Script name - `regrid_popden_files.py`
* Create multiple population density files for E3SM-GCAM at different resolutions
* One input argument determines resolution