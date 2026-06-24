# Scripts for generating baseline data
--------------------------------------------------

## Create baseline `npp`, `hr`, and `pft_wt` files for E3SM-GCAM simulations:
### Script name - `create_e3sm_gcam_land_scalar_baseline_local.r`

* This function reads monthly elm history files to obtain pft-level values for baseline of:
	* `veg_cf%npp`
	* `col_cf%hr`
	* `veg_pp%wtgcell` for the veg landunit only
	* total cell area
* This data is currently stored in the `h2` file in our current run script.

--------------------------------------------------

## Create baseline `hdd`, and `cdd` files for E3SM-GCAM simulations:
### Script name - `era5_hdd_cdd_elm_grids.py`

 Computes mean annual Heating Degree Days (HDD) and Cooling Degree Days (CDD)
 from ERA5 near-surface temperature forcing, then regrids the spatial fields
 to each of the five standard ELM grid resolutions used in E3SM-GCAM coupling.
