# Author - Dalei Hao, Pacific Northwest National Lab

# iesm2_land_mapping.r

# generate new glm input files and iesmv2 iac grid mapping files that reflect current gcam regions/glus
# this is for iesm in e3sm v2, using the original GLM

# the glm grid files are at half-degree, and only one land unit per grid cell can be assigned
# so use the dominant land unit per grid cell

# in order to be compatible with GLM, the region and glu numbers need to be an enumeration up to the region number
#	and maximum glu number in each region
# this is based on the available gcam xml output list, not the moirai regionXglu combos

# the moirai data are from v3.0, which is used as the land basis for iESM v2 with GCAM 5.1

# all gridded data have first value at upper left corner -180,90
# except that netcdf files have their own formats

# GCAM now outputs Taiwan_Taiwan and Taiwan_ChinaCst
# so there are 32 regions, and the taiwan cells do not have to be remapped to china

# required libraries

library(raster)
library(rasterVis)
#library(ggplot2)
#library(maptools)
#library(ggmap)
#library(ggplot2)
#library(maps)
library(rgdal)
library(ncdf4)
library(XML)
library(parallel)
library(plyr)

#####
# function: iesm2_land_mapping(input_dir, new_dir, write_elm2gcam = TRUE, write_gcam2glm = TRUE)
#
# four arguments
# input_dir:	the directory containing the 15 input files (not counting headers and other aux files), see within function
# new_dir:		the directory to write the 14 output files to
# write_elm2gcam:	TRUE = generate the elm2gcam grid mapping file (takes about 4 hours on my desktop)

# 14 output files:

# new region map, codes, names, and continents
# the region codes, names, and continent codes are in matching order
# region and zone codes need to be enumerated based on the gcam xml order

# gcam_region_grid.txt:				half-degree gcam region map as a text file with no header, codes enumerated by woodharvest.xml order
# gcam_region_grid.asc:				half-degree gcam region map as an esri ascii grid - not used by iesm, but useful for viewing
# gcam_region_codes.txt:			gcam region codes, as enumerated by woodharvest.xlm order; used in gcam_region_grid
# gcam_region_names.txt:			gcam region names, in order to match the codes in gcam_region_codes.txt
# continent_codes_region.txt:		glm continent code for each gcam region, in the order of gcam_region_codes.txt
# gcam_zone_grid.txt:				half-degree gcam glu map as a text file with no header, codes enumerated within each region by woodharvest.xml order
# gcam_zone_grid.txt:				half-degree gcam glu map as an esri ascii grid
# reg2ctry_mapping.txt:				gcam region code and corresponding glm country code for eahc glm country, in order of cnames.txt.sort2wh
# new_continent.codes.txt.sort2wh:	glm continent code for each glm country, in order of cnames.txt.sort2wh
# vba_LUH1format.nc:				netcdf glm biomass file, generated from vba_LUH1format.txt
# initial_state_LUH2_2015_v3.nc:	glm initial state updated with glm cell area, updated from initial_state_LUH2_2015_v2.nc
# gcam2glm_mapping.csv:				grid mapping of gcam regionXglu to glm half-degree grid, with weights
# elm0.9x1.25togcam_mapping.csv:	grid mapping of gcam regionXglu to elm nominal 1-degree grid, with weights
# iac_region_glu_codes.csv:			diagnostic; out region and glu names and codes for iac, as per the order in woodharvest.xml
#										out_reg_code and out_basin_code are used in gcam_region_grid and gcam_zone_grid for glm; out_lu_code is used in the iac


######### iesm2_land_mapping function definition
# two helper functions for the grid mappings are defined below

input_dir = "/Users/haod776/Library/CloudStorage/OneDrive-PNNL/Documents/work/E3SM/GCAM/mapping/input"
new_dir = "/Users/haod776/Library/CloudStorage/OneDrive-PNNL/Documents/work/E3SM/GCAM/mapping/output"
write_elm2gcam = TRUE

cat("Starting iesm2_land_mapping", date(), "\n")

# make sure there is a "/" at the end of input path
pl = nchar(input_dir)
if(substr(input_dir,pl,pl) != "/") {
	input_dir = paste0(input_dir,"/")
}
	
# make sure there is a "/" at the end of output path
plo = nchar(new_dir)
if(substr(new_dir,plo,plo) != "/") {
	new_dir = paste0(new_dir,"/")
}

# create output directory
dir.create(new_dir)

num_cores = detectCores()

# projection info for all files, data, and rasters
PROJ4_STRING = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# and extent object for setting raster object extent to the full globe
FULLGLOBE = extent(-180,180,-90,90)


#china_code = 11			# this is the moirai country code
#china_name = "China"
#taiwan_code = 30			# this is the moirai county code
#taiwan_name = "Taiwan" # this is the same name for the region and the basin
#taiwan_basin_code = 103
#china_coast_basin_code = 78 # this is the other basin in taiwan, only 12 cells

# new elm to gcam mapping file
elm2gcam_on = paste0(new_dir, "elm0.9x1.25tocountry_mapping.csv")


##### now read in all relevant input data

## current gcam regionXglu raster map name
# the values are country code
#moirai_fn = paste0(input_dir, "country_out.bil")
#moirai_rast <- raster(moirai_fn)

moirai_fn = paste0(input_dir, "country_out.bil")
moirai_in <- readBin(moirai_fn, what="numeric", n=2160*4320, size = 4)
ndinds = which(moirai_in == -9999)
moirai_in[ndinds] = NA
dim(moirai_in) <- c(4320,2160)
moirai = t(moirai_in)
moirai_rast = raster(x= moirai,xmn=-180,ymn=-90,xmx = 180,ymx=90, crs=PROJ4_STRING)


## current gcam cell area raster map name
# this corresponds with the moirai regionXglu map
# but need to extrapolate cell area outside of land mask to get full coarse cell area
# so fill NA cell areas by latitide
# this takes about 20 minutes
moirai_cell_area_fn = paste0(input_dir, "hyde_cell_plus.bil")
moirai_cell_area_in <- readBin(moirai_cell_area_fn, what="numeric", n=2160*4320, size = 4)
ndinds = which(moirai_cell_area_in == -9999)
moirai_cell_area_in[ndinds] = NA
dim(moirai_cell_area_in) <- c(4320,2160)
moirai_cell_area = t(moirai_cell_area_in)
moirai_cell_area_rast = raster(x= moirai_cell_area,xmn=-180,ymn=-90,xmx = 180,ymx=90, crs=PROJ4_STRING)
for(i in 1: 2160){
	nainds = is.na(moirai_cell_area_rast[i,])
	moirai_cell_area_rast[i,][nainds] = min(moirai_cell_area_rast[i,], na.rm=TRUE)
}

## current gcam land area raster map name
# this corresponds with the moirai regionXglu map
moirai_land_area_fn = paste0(input_dir, "hyde_land_plus.bil")
moirai_land_area_in <- readBin(moirai_land_area_fn, what="numeric", n=2160*4320, size = 4)
ndinds = which(moirai_land_area_in == -9999)
moirai_land_area_in[ndinds] = NA
dim(moirai_land_area_in) <- c(4320,2160)
moirai_land_area = t(moirai_land_area_in)
moirai_land_area_rast = raster(x= moirai_land_area,xmn=-180,ymn=-90,xmx = 180,ymx=90, crs=PROJ4_STRING)


## original clm .9x1.25 surface file for grid specification
# need LONE, LONW, LATN, LATS, by (lsmlat,lsmlon), and probably lsmlat and lsmlon
# read in below for elm2gcam mapping
clm_surf_fn = paste0(input_dir, "surfdata_0.9x1.25_ZGICN32c_c120807.nc")


##### pre-first
# the gcam output list determines the region and zone values by order
# so build a table with all the relevant info
# the gcam codes are the moirai values
# the out codes are based on the gcam output order
# gcam output must be grouped by region

# gcam moirai region list
country_fn = paste0(input_dir, "FAO_iso_VMAP0_ctry.csv")
country_in = read.csv(country_fn, sep=",", header=TRUE);

gcam_lu = country_in
colnames(gcam_lu) = c('fao_code','iso3_abbr','fao_name','gcam_cty_code','gcam_cty_name')



############# elm2gcam single-cell grid mapping function

# needs cell index and four arguments: moirai_rast, moirai_cell_area_rast, moirai_land_area_rast, gcam_lu
elm2gcam_cell <- function(cell, regglu_rast, cell_area_rast, land_area_rast, gcam_lu_df)	{
  # this assumes that cell starts at 0 to eliminate a calc inside
  lt = trunc(cell / num_lon) + 1
  ln = cell %% num_lon + 1
  
  lnmin = lonW[ln,lt]
  lnmax = lonE[ln,lt]
  ltmin = latS[lt,lt]
  ltmax = latN[lt,lt]
  
  # need to shift the longitude values so that they range from -180 to 180 rather than -0.625 to 179.375
  # deal with the boundary straddle also
  if (lnmin < 180 & lnmax > 180) {
    # need two spatial objects
    # need each conrer of the larger cell, clockwise, and first and last are the same
    lnmax_tmp = lnmax
    lnmax = 180
    corners = cbind(c(lnmin, lnmax, lnmax, lnmin, lnmin), c(ltmax, ltmax, ltmin, ltmin, ltmax))
    lnmin = -180
    lnmax =  -180 + (lnmax_tmp - 180)
    corners2 = cbind(c(lnmin, lnmax, lnmax, lnmin, lnmin), c(ltmax, ltmax, ltmin, ltmin, ltmax))
    p1 = Polygons( list(Polygon(coords = corners, hole=FALSE)) , 1)
    p2 = Polygons( list(Polygon(coords = corners2, hole=FALSE)) , 2)
    spobj = SpatialPolygons( list(p1, p2) , proj4string=CRS(PROJ4_STRING))
  } else { # only one spatial object needed
    if (lnmin >= 180 & lnmax > 180) {
      lnmin = -180 + (lnmin - 180)
      lnmax = -180 + (lnmax - 180)
    }
    # need each conrer of the larger cell, clockwise, and first and last are the same
    corners = cbind(c(lnmin, lnmax, lnmax, lnmin, lnmin), c(ltmax, ltmax, ltmin, ltmin, ltmax))
    spobj = SpatialPolygons( list(Polygons( list(Polygon(coords = corners, hole=FALSE)) , 1)) , proj4string=CRS(PROJ4_STRING))
  }
  
  
  lu_cells = extract(regglu_rast, spobj, weights=TRUE, normalizeWeights=FALSE, cellnumbers=TRUE, df=TRUE, small=TRUE, exact=TRUE)
  names(lu_cells)[names(lu_cells) == "layer"] = "gcam_cty_code"
  ca = extract(cell_area_rast, spobj, weights=TRUE, normalizeWeights=FALSE, cellnumbers=TRUE, df=TRUE, small=TRUE, exact=TRUE)
  la = extract(land_area_rast, spobj, weights=TRUE, normalizeWeights=FALSE, cellnumbers=TRUE, df=TRUE, small=TRUE, exact=TRUE)
  names(ca)[names(ca) == "layer"] = "cell_area"
  names(la)[names(la) == "layer"] = "land_area"
  lu_cells=merge(lu_cells, ca, by=c("ID", "cell"))
  lu_cells=merge(lu_cells, la, by=c("ID","cell"))
  lu_cells$weighted_land_area = lu_cells$weight * lu_cells$land_area
  lu_cells$weighted_cell_area = lu_cells$weight * lu_cells$cell_area
  # need area of entire coarse cell, but fine cells have NA where there is no land
  total_cell_area = sum(lu_cells$weighted_cell_area, na.rm=TRUE)
  # filter out empty values and zero area
  lu_cells = lu_cells[!is.na(lu_cells$gcam_cty_code) & lu_cells$weighted_land_area > 0 & !is.na(lu_cells$weighted_land_area),]
  total_land_area = sum(lu_cells$weighted_land_area, na.rm=TRUE)
  
  
  if (total_land_area > 0 & total_cell_area > 0) {
    
    lu_area = aggregate(weighted_land_area ~ gcam_cty_code, lu_cells, FUN=sum, na.rm=TRUE)
    lu_area$total_land_area = total_land_area
    lu_area$total_cell_area = total_cell_area
    lu_area$lu_fraction = lu_area$weighted_land_area / lu_area$total_cell_area
    #lu_area$gcam_lu_code[lu_area$gcam_lu_code == taiwan_code * 10000 + taiwan_basin_code | gcam_lu_code == taiwan_code * 10000 + china_coast_basin_code] =
    #	gcam_lu_code = china_code * 10000 + taiwan_code
    weight = lu_area$lu_fraction
    gcam_cty_code = lu_area$gcam_cty_code
    numrec = length(weight)
    xlon = rep(ln, numrec)
    ylat = rep(lt, numrec)
    df = data.frame(gcam_cty_code, xlon, ylat, weight)
    df = merge(df, gcam_lu_df, by=c("gcam_cty_code"), sort=FALSE)
    df_out = df[,c("gcam_cty_code","gcam_cty_name", "xlon", "ylat", "weight")]
    return(df_out)
  } # end if values exist in this cell
} # end single sell elm2gcam map function


# do this only if asked
if (write_elm2gcam) {

	clm_surf = nc_open(clm_surf_fn)
	num_lon = clm_surf$dim$lsmlon$len
	num_lat = clm_surf$dim$lsmlat$len
	#lsmlon = ncvar_get(clm_surf,varid="lsmlon",start=c(1), count=c(num_lon))
	#lsmlat = ncvar_get(clm_surf,varid="lsmlat",start=c(1), count=c(num_lat))
	lonE = ncvar_get(clm_surf,varid="LONE",start=c(1,1), count=c(num_lon, num_lat))
	lonW = ncvar_get(clm_surf,varid="LONW",start=c(1,1), count=c(num_lon, num_lat))
	latN = ncvar_get(clm_surf,varid="LATN",start=c(1,1), count=c(num_lon, num_lat))
	latS = ncvar_get(clm_surf,varid="LATS",start=c(1,1), count=c(num_lon, num_lat))
	nc_close(clm_surf)

	num_cells = num_lon * num_lat

	# run the single-cell elm2gcam grid mapping function in parallel - this takes about 50 minutes
	cat("\nstarting elm2gcam mapping", date(), "\n")
	mcout = mclapply(c(0:(num_cells-1)), function(i) elm2gcam_cell(cell = i, moirai_rast, moirai_cell_area_rast, moirai_land_area_rast, gcam_lu), mc.cores = num_cores)
	cat("\nfinishing elm2gcam2 mapping (before sorting and writing)", date(), "\n")

	# need to merge the data frames
	elm2gcam_out = rbind.fill(mcout)

	# now sort the table by land unit and location and write it
	elm2gcam_out = elm2gcam_out[order(elm2gcam_out$gcam_cty_name, elm2gcam_out$xlon, elm2gcam_out$ylat),]
	
	# change the column names
	names(elm2gcam_out) <- c("country_index", "gcam_cty_name", "x", "y", "Weight")
	
	write.table(elm2gcam_out, elm2gcam_on, row.names=FALSE, quote=FALSE, sep=",")

} # end if write_elm2gcam = TRUE


