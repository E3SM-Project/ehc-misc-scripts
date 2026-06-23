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

