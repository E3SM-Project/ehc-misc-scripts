# Author - Alan Di Vittorio, Lawrence Berkeley National Lab

# e3sm_gcam_land_mapping.r

# running the function does not allow the parallel functions to work
#    run this as a script by directly running the code within the function

# generate new glm input files and iesmv2 iac grid mapping files that reflect current gcam regions/glus
# this is for iesm in e3sm v2 and v3, using the original GLM
# now both the region elm2gcam and usa-state elm2gcam files are written

# the glm grid files are at half-degree, and only one land unit per grid cell can be assigned for id thematics
# so use the dominant land unit per grid cell

# the two mapping files are fraction of output cell, except for a bug that does not include the full cell area for normailization

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
#library(rgdal)
library(ncdf4)
#library(XML)
library(xml2)
library(parallel)
library(plyr)
library(terra)

#####
# function: e3sm_gcam_land_mapping.r(input_dir, new_dir, write_elm2gcam = TRUE, write_gcam2glm = TRUE)
#
# four arguments
# input_dir:	the directory containing the 15 input files (not counting headers and other aux files), see within function
# new_dir:		the directory to write the 14 output files to
# write_elm2gcam:	TRUE = generate the elm2gcam grid mapping file (takes about 30 minutes for three files on my desktop for f09)
# write_gcam2glm	:   TRUE = generate the gcam2glm grid mapping file (takes about 25 minutes on my desktop for f09)
# the whole thing took about 55 minutes for f09

# 26 output files (unless write_elm2gcam or write_gcam2glm are FALSE below):

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
# iac_region_glu_codes.csv:			diagnostic; out region and glu names and codes for iac, as per the order in woodharvest.xml
#										out_reg_code and out_basin_code are used in gcam_region_grid and gcam_zone_grid for glm;
#                                       out_lu_code is used in the iac
#
# elm0.9x1.25togcam_mapping.csv:	grid mapping of gcam regionXglu to elm nominal 1-degree grid, with weights
# elm1.9x2.5togcam_mapping.csv:	grid mapping of gcam regionXglu to elm nominal 2-degree grid, with weights
# elm0.5x0.5togcam_mapping.csv:	grid mapping of gcam regionXglu to elm nominal half-degree grid, with weights
# elm0.125x0.125togcam_mapping.csv:	grid mapping of gcam regionXglu to elm nominal eighth-degree grid, with weights
# elm0.25x0.25togcam_mapping.csv:	grid mapping of gcam regionXglu to elm nominal quarter-degree grid, with weights
# note that these files map cell indices in the order and origin defined by each grid (values increase, lon varies fastest)
# grids < 1-deg have cells aligned with -180, -90 as the edge of the first cell
# grids >= 1-deg have cells aligned with 0, -90 as the center of the first cell

# also five elm##x##togcam_usa_mapping.csv files, for the five resolutions
#    and iac_state_glu_codes.csv (diagnostic), gcam_state_region_codes.txt, gcam_state_region_names.txt  

# should be able to just add another input and output file below to get another resolution elm2gcam mapping file

# as noted above, the parallel function have not been working when calling the function,
#    so set the arguments here and run the code in the function by lines
#    need to set inputdata_dir to the local inputdata directory
inputdata_dir = "/global/cfs/cdirs/e3sm/inputdata"
input_dir = paste0(inputdata_dir, "/iac/giac/inputs_for_mapping_generation")
new_dir = "./mapping_outfiles"
write_elm2gcam = TRUE
write_gcam2glm = FALSE

######### e3sm_gcam_land_mapping function definition
# two helper functions for the grid mappings are defined below

e3sm_gcam_land_mapping <- function(input_dir = paste0(inputdata_dir, "/iac/giac/inputs_for_mapping_generation"),
		new_dir = "./mapping_outfiles", write_elm2gcam = TRUE, write_gcam2glm = FALSE) {

cat("Starting e3sm_gcam_land_mapping", date(), "\n")

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

GLM_NODATA = 0

# this is to check for precision error when calculating the elm cell boundaries
BTOL = 1e-6

#china_code = 11			# this is the moirai country code
#china_name = "China"
#taiwan_code = 30			# this is the moirai county code
#taiwan_name = "Taiwan" # this is the same name for the region and the basin
#taiwan_basin_code = 103
#china_coast_basin_code = 78 # this is the other basin in taiwan, only 12 cells

# set the GCAM output states and their abbreviations (region names)
gcam_states_names = c("Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado", "Connecticut", "Delaware", "District of Columbia",
	"Florida", "Georgia", "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa", "Kansas", "Kentucky", "Louisiana", "Maine", "Maryland", "Massachusetts",
	"Michigan", "Minnesota", "Mississippi", "Missouri", "Montana", "Nebraska", "Nevada", "New Hampshire", "New Jersey", "New Mexico", "New York",
	"North Carolina", "North Dakota", "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota", "Tennessee",
	"Texas", "Utah", "Vermont", "Virginia", "Washington", "West Virginia", "Wisconsin", "Wyoming")
gcam_states_abr = c("AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE", "DC",
	"FL", "GA", "HI", "ID", "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD", "MA",
	"MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY",
	"NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN",
	"TX", "UT", "VT", "VA", "WA", "WV", "WI", "WY")
gcam_states_df = data.frame(gcam_states_names, gcam_states_abr)

## new file names

# new region map, codes, names, and continents
# the region codes, names, and continent codes are in matching order
# region and zone codes need to be enumerated based on the gcam xml order
reg_on = paste0(new_dir, "gcam_region_grid.txt")
reg_an = paste0(new_dir, "gcam_region_grid.asc")
regcode_on = paste0(new_dir, "gcam_region_codes.txt")
regname_on = paste0(new_dir, "gcam_region_names.txt")
contcode_reg_on = paste0(new_dir, "continent_codes_region.txt")
state_regcode_on = paste0(new_dir, "gcam_state_region_codes.txt")
state_regname_on = paste0(new_dir, "gcam_state_region_names.txt")

# new zone map, this is based on the maximum number of gcam output glus in a region, and enumerated
# the actual glu numbers are not used
zone_on = paste0(new_dir, "gcam_zone_grid.txt")
zone_an = paste0(new_dir, "gcam_zone_grid.asc")

# new region to glm country mapping
reg2ctry_on = paste0(new_dir, "reg2ctry_mapping.txt")

# new continent codes sorted to glm countries
contcode_ctry_on = paste0(new_dir, "new_continent.codes.txt.sort2wh")

# new gcam to glm luc mapping file (for gcam2glm)
gcam2glm_on = paste0(new_dir, "gcam2glm_mapping.csv")

# new netcdf glm biomass file name (for gcam2glm)
bio_on = paste0(new_dir, "vba_LUH1format.nc")

######## new elm to gcam mapping files

# these are the standard grid mapping files used for the land
elm2gcam1_on = paste0(new_dir, "elm0.9x1.25togcam_mapping.csv")
elm2gcam2_on = paste0(new_dir, "elm1.9x2.5togcam_mapping.csv")
elm2gcamh_on = paste0(new_dir, "elm0.5x0.5togcam_mapping.csv")
elm2gcame_on = paste0(new_dir, "elm0.125x0.125togcam_mapping.csv")
elm2gcamq_on = paste0(new_dir, "elm0.25x0.25togcam_mapping.csv")
elm2gcam_ons = c(elm2gcam1_on, elm2gcam2_on, elm2gcamh_on, elm2gcame_on, elm2gcamq_on)

# these are a new set, based on the land ones that are used for co2 mapping and include states
elm2gcam1_usa_on = paste0(new_dir, "elm0.9x1.25togcam_usa_mapping.csv")
elm2gcam2_usa_on = paste0(new_dir, "elm1.9x2.5togcam_usa_mapping.csv")
elm2gcamh_usa_on = paste0(new_dir, "elm0.5x0.5togcam_usa_mapping.csv")
elm2gcame_usa_on = paste0(new_dir, "elm0.125x0.125togcam_usa_mapping.csv")
elm2gcamq_usa_on = paste0(new_dir, "elm0.25x0.25togcam_usa_mapping.csv")
elm2gcam_usa_ons = c(elm2gcam1_usa_on, elm2gcam2_usa_on, elm2gcamh_usa_on, elm2gcame_usa_on, elm2gcamq_usa_on)


# new current 2015 initial glm file, with cell_area added for gcam2glm
init_glm_on = paste0(new_dir, "initial_state_LUH2_2015_v3.nc")

iac_rg_codes_on = paste0(new_dir, "iac_region_glu_codes.csv")
iac_state_codes_on = paste0(new_dir, "iac_state_glu_codes.csv")

##### now read in all relevant input data

## current gcam regionXglu raster map name
# the values are region code * 10000 + glu
moirai_fn = paste0(input_dir, "moirai_valid_region32_water_basin235.bsq")
moirai_rast <- raster(moirai_fn)


## current gcam state delineation
# this is just for co2 emissions grid mapping
# create the same format as the land grid mapping, even though all the info isn't needed, so that the read function has minimal changes
# add 50 to the state ids so as to not conflict with the region ids
# the number is not important in the output, just the names, which are the abbreviations
# these will be new region ids for the us states
# need to create the raster codes below, which are: region code * 10000 + glu
# and add the unique records to the gcam land_unit list of regionXglus

gcam_usa_state_fn = paste0(input_dir, "gcamusa_state_glu_wgs84/gcamusa_state_glu_wgs84.shp")
state_vector = terra::vect(gcam_usa_state_fn)
moirai_terra <- terra::rast(moirai_fn)
in_state_terra = terra::rasterize(state_vector, moirai_terra, field="state_id", fun=min, update=FALSE, background=NA)
in_state_rast = raster(in_state_terra)
in_state_rast = deratify(in_state_rast, att="state_id", fun="min")
in_state_rast = in_state_rast + 50


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

## current gcam regionXglu combinations
# this contains the region names (with spaces) and their codes, with glu codes (newer versions will include glu names, but these are not used here)
moirai_codes_fn = paste0(input_dir, "MOIRAI_reggcam_GLU_with_spaces.csv")
moirai_codes = read.csv(moirai_codes_fn, skip=4, header=TRUE, stringsAsFactors=FALSE)

## current gcam basin (glu) to country mapping with basin abbreviations for woodharvest.xml
# only the GCAM_basin_ID (just the number), GLU_name (the abbreviation) and the Basin_name (the full name) are needed
gcam_basin2ctry_fn = paste0(input_dir, "basin_to_country_mapping.csv")
gcam_basin2ctry = read.csv(gcam_basin2ctry_fn, skip=7, header=TRUE, stringsAsFactors=FALSE)

## current xml woodharvest file gcam land unit = regionXbasin in iesm
# this has same list as luc.xml, but simpler to parse
# these are the current regionXglu combos in iesm (names are regionXglu(abr.) combined, and spaces removed from region names)
# there are currently 392, even though only 384 are really valid in GCAM
# because every combo does not have gcam output
# these appear to be the same land units as for luc
# remove whitespace for processing
# create the data frame below
wh_xml_fn = paste0(input_dir, "woodharvest.xml")
wh_xml = read_xml(readChar(wh_xml_fn, file.info(wh_xml_fn)$size))
# the land units are the second child, and the first record is a comment
# xml_length does not count the comment record
wh_lu_child = xml_child(wh_xml,search=2)
num_lus = xml_length(wh_lu_child)
land_unit=array(dim=num_lus)
wh_lu_contents = xml_contents(wh_lu_child)
for(i in 2:(num_lus+1)){
	land_unit[i-1] = xml_text(wh_lu_contents[i])
}


#wh_xml = xmlTreeParse(wh_xml_fn, useInternal=TRUE)
#top = xmlRoot(wh_xml)
#text=xmlSApply(top,function(x) xmlSApply(x,xmlValue)) # want text$column for land units, and the first element is just a comment
#land_unit = unlist(text$column[2:length(text$column)])


# current 2015 glm text biomass file (kgC/m^2)
bio_fn = paste0(input_dir, "vba_LUH1format.txt")
bio_in = read.csv(bio_fn, sep=" ", dec=".", header=FALSE)
bio_mat = data.matrix(bio_in)
bio_rast = raster(x= bio_mat,xmn=-180,ymn=-90,xmx = 180,ymx=90, crs=PROJ4_STRING)

# current 2015 initial glm file
# read in below for gcam2glm mapping
init_glm_fn = paste0(input_dir, "initial_state_LUH2_2015_v2.nc")

## elm surfdata files for grid specification
## replace this with more current files; must use LONGXY and LATIXY
## original clm .9x1.25 surface file for grid specification
# need LONE, LONW, LATN, LATS, by (lsmlat,lsmlon), and probably lsmlat and lsmlon
# read in below for elm2gcam mapping
#clm_surf_fn = paste0(input_dir, "surfdata_0.9x1.25_ZGICN32c_c120807.nc")
clm_surf1_fn = paste0(input_dir, "surfdata_0.9x1.25_HIST_simyr2015_c201021.nc")
clm_surf2_fn = paste0(input_dir, "surfdata_1.9x2.5_SSP5_RCP85_simyr2015_c210916.nc")
clm_surfh_fn = paste0(input_dir, "surfdata_0.5x0.5_HIST_simyr2015_c220318.nc")
clm_surfe_fn = paste0(input_dir, "surfdata_0.125x0.125_HIST_simyr2015_c241205.nc")
clm_surfq_fn = paste0(input_dir, "surfdata_0.25x0.25_simyr2015_c250312.nc")
clm_surf_fns = c(clm_surf1_fn, clm_surf2_fn, clm_surfh_fn, clm_surfe_fn, clm_surfq_fn)


## original file names for inputs

# orginal other initialization file - to get cell area
orig_other_fn = paste0(input_dir, "gothr_1500-2005.nc")

# original country map half degree
ctry_fn = paste0(input_dir, "ccodes_half_deg.txt")
ctry_in = read.csv(ctry_fn, sep=" ", header=FALSE)
ctry_rast = data.matrix(ctry_in)
ctry_rast = raster(x= ctry_rast,xmn=-180,ymn=-90,xmx = 180,ymx=90, crs=PROJ4_STRING)

# original ctry codes, sorted to wood harvest
ctrycode_fn = paste0(input_dir, "ccodes.txt.sort2wh")
ctrycode = read.csv(ctrycode_fn, header=FALSE, stringsAsFactors=FALSE)

# original ctry names, sorted to wood harvest
ctryname_fn = paste0(input_dir, "cnames.txt.sort2wh")
ctryname = read.csv(ctryname_fn, header=FALSE, stringsAsFactors=FALSE, sep="")

# original continent map half degree
cont_fn = paste0(input_dir, "gcodes_continent_half_deg_DUMMY.asc")
cont_in = read.csv(cont_fn, sep=" ", header=FALSE)
cont_rast = data.matrix(cont_in)
cont_rast = raster(x= cont_rast,xmn=-180,ymn=-90,xmx = 180,ymx=90, crs=PROJ4_STRING)

# original netcdf biomass file for gcam2glm
# variable is float biomass[lon,lat,time] (R order); time is size 1
old_bio_fn = paste0(input_dir, "miami_biomass_conform_0.5x0.5_map.nc")
old_bio = nc_open(old_bio_fn)
nbio_lon = old_bio$dim$lon$len
nbio_lat = old_bio$dim$lat$len
nc_close(old_bio)

# the following files are original glm inputs that are not needed here
# the read in code remains in case comparison is needed

# path to original glm inputs
#orig_dir = "~/projects/e3sm/giac_v2/glm_inputs/original/"

# original region map half degree - stored lat-line by line
#reg_fn = paste0(orig_dir, "AEZ_region_grid.txt")
# raster assumes that the upper left corner is stored in matrix[1,1], and that matrix rows correspond with latitude
#reg_in = read.csv(reg_fn, sep=" ", dec=".", header=FALSE)
#reg_rast = data.matrix(reg_in)
#reg_rast = raster(x=reg_rast,xmn=-180,ymn=-90,xmx = 180,ymx=90, crs=PROJ4_STRING)

# original gcam region codes
#reg14code_fn = paste0(orig_dir, "codes_halfdeg_minicam.txt")
#reg14code = read.csv(reg14code_fn, header=FALSE, stringsAsFactors=FALSE)

# original gcam region names, matching the region codes
#reg14name_fn = paste0(orig_dir, "names_minicam.txt")
#reg14name = read.csv(reg14name_fn, header=FALSE, stringsAsFactors=FALSE)

# original continent codes matching original gcam region codes/names
#contcode_reg_fn = paste0(orig_dir, "continent_codes_minicam_test.txt")
#contcode_reg = read.csv(contcode_reg_fn, header=FALSE, stringsAsFactors=FALSE)

# original gcam region to glm country mapping (turkey is in ee)
#reg2ctry_fn = paste0(orig_dir, "codes2glm_gcam_turkey_in_ee.txt")
#reg2ctry = read.csv(reg2ctry_fn, header=FALSE, stringsAsFactors=FALSE)

# original aez map half degree
#aez_fn = paste0(orig_dir, "AEZ_zone_grid.txt")
#aez_in = read.csv(aez_fn, sep=" ", dec=".", header=FALSE)
#aez_rast = data.matrix(aez_in)
#aez_rast = raster(x=aez_rast,xmn=-180,ymn=-90,xmx = 180,ymx=90, crs=PROJ4_STRING)

# original continent codes, sorted to wood harvest country codes/names
#contcode_ctry_fn = paste0(orig_dir, "continent.codes.txt.sort2wh")
#contcode_ctry = read.csv(contcode_ctry_fn, header=FALSE, stringsAsFactors=FALSE)

# this is not used in this version of GLM, so don't read it in or worry about it
# this is the same as the AEZ_region_grid.txt file above
#greg_fn = paste0(orig_dir, "regcodes_halfdeg.txt")
#greg_in = read.csv(greg_fn, sep=" ", header=FALSE)
#greg_rast = data.matrix(greg_in)
#greg_rast = raster(x= greg_rast,xmn=-180,ymn=-90,xmx = 180,ymx=90, crs=PROJ4_STRING)


##### pre-first
# the gcam output list determines the region and zone values by order
# so build a table with all the relevant info
# the gcam codes are the moirai values
# the out codes are based on the gcam output order
# gcam output must be grouped by region

# gcam moirai region list
reg_list = unique(moirai_codes[,c("gcam_reg_code","gcam_reg_name")])
num_reg = nrow(reg_list)

# create df and separate region and glu in the gcam output xml list
gcam_lu = data.frame(land_unit)
gcam_lu$land_unit = unlist(as.character(gcam_lu$land_unit))
gcam_lu$out_lu_code = c(1:nrow(gcam_lu)) # enumerated land units by order
# region spaces are now in the input file so don't remove them
#gcam_lu$land_unit = unlist(lapply(gcam_lu$land_unit, function(x) gsub("\\s+","",x)))
gcam_lu$gcam_reg_name = NA
gcam_lu$gcam_reg_code = NA
gcam_lu$basin_abr = NA
for(r in 1:nrow(reg_list)){
	gcam_lu$gcam_reg_name[grep(paste0("^",reg_list$gcam_reg_name[r]),gcam_lu$land_unit)] = reg_list$gcam_reg_name[r]
	gcam_lu$gcam_reg_code[grep(paste0("^",reg_list$gcam_reg_name[r]),gcam_lu$land_unit)] = reg_list$gcam_reg_code[r]
	gcam_lu$basin_abr[grep(paste0("^",reg_list$gcam_reg_name[r]),gcam_lu$land_unit)] = 
		substr(sapply(regmatches(gcam_lu$land_unit[grep(paste0("^",reg_list$gcam_reg_name[r]),gcam_lu$land_unit)], 
			regexpr(paste0("^",reg_list$gcam_reg_name[r]),gcam_lu$land_unit[grep(paste0("^",reg_list$gcam_reg_name[r]),gcam_lu$land_unit)]), invert=TRUE), "[[", 2),
				2,
				nchar(sapply(regmatches(gcam_lu$land_unit[grep(paste0("^",reg_list$gcam_reg_name[r]),gcam_lu$land_unit)],
					regexpr(paste0("^",reg_list$gcam_reg_name[r]),gcam_lu$land_unit[grep(paste0("^",reg_list$gcam_reg_name[r]),gcam_lu$land_unit)]), invert=TRUE), "[[",2)))
}

# compile table of unique out regions and write the region code and region name files
# the moirai and out gcam regions should be the same
gcam_reg_name = unique(gcam_lu[,c("gcam_reg_name")])
out_reg_list = data.frame(gcam_reg_name, stringsAsFactors=FALSE)
num_out_reg = nrow(out_reg_list)
out_reg_list$out_reg_code = c(1:num_out_reg)
write.table(out_reg_list$out_reg_code, regcode_on, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(out_reg_list$gcam_reg_name, regname_on, row.names=FALSE, col.names=FALSE, quote=FALSE)

# add the out region codes to the lu table
gcam_lu = merge(gcam_lu, out_reg_list, by = c("gcam_reg_name"), sort=FALSE)

# merge the full basin names and gcam basin codes here
gcam_lu = merge(gcam_lu, gcam_basin2ctry[,c("GCAM_basin_ID", "Basin_name", "GLU_name")], by.x = c("basin_abr"), by.y = c("GLU_name"), sort=FALSE)
names(gcam_lu)[which(names(gcam_lu) == "GCAM_basin_ID")] = "gcam_basin_code"

gcam_lu = gcam_lu[order(gcam_lu$out_lu_code),]

# loop through the regions to enumerate the basins; store and print the max basins per region
gcam_lu$out_basin_code = NA
max_basins = 0
mb_reg_ind = NA
out_reg_list$num_basin = NA
for (r in 1:num_out_reg) {
	num_temp_basin = length(gcam_lu$out_basin_code[gcam_lu$out_reg_code == r])
	gcam_lu$out_basin_code[gcam_lu$out_reg_code == r] = c(1:num_temp_basin)
	out_reg_list$num_basin[out_reg_list$out_reg_code == r] = num_temp_basin
	if (num_temp_basin > max_basins) {
		max_basins = num_temp_basin
		mb_reg_ind = r
	}
}
cat("\nMax output basins, within region", mb_reg_ind, out_reg_list$gcam_reg_name[mb_reg_ind], "is", max_basins, "\n")

# write the diagnostic output gcam region and glu codes and names

iac_gcam_reg_glu_codes = gcam_lu[,c("out_lu_code", "land_unit","out_reg_code", "out_basin_code", "gcam_reg_name", "Basin_name",
									"basin_abr", "gcam_reg_code", "gcam_basin_code")]
iac_gcam_reg_glu_codes = iac_gcam_reg_glu_codes[order(iac_gcam_reg_glu_codes$out_lu_code),]
write.table(iac_gcam_reg_glu_codes, iac_rg_codes_on, row.names=FALSE, quote=FALSE, sep=",")


#####################################
######## create different region data for gcam-usa grid mapping
# this will include overlapping USA and state regions; the will be parsed by the E3SM-GCAM code

## add the states to gcam_lu
# add 50 to the state ids to ensure they are different that the region ids (the raster has already been adjusted in this way)
# need to create all of the gcam_lu columns
# gcam-usa 51 state region names are the abbreviations, plus DC
usa_df = terra::as.data.frame(state_vector)
usa_df$key = NULL
usa_df = merge(usa_df, unique(gcam_lu[,c("basin_abr", "gcam_basin_code")]), by.x = "glu_id", by.y = "gcam_basin_code", all.x = TRUE)
usa_df$state_id = as.integer(usa_df$state_id) + 50
usa_df$land_unit = paste0(usa_df$state_nm, "_", usa_df$basin_abr)
usa_df$out_lu_code = NA
usa_df$out_basin_code = NA

# state names are abbrviations; 50 states plus DC
state_df = unique(usa_df[,c("state_nm", "state_id")])
state_df = merge(state_df, gcam_states_df, by.x = c("state_nm"), by.y = c("gcam_states_names"))
# the merge automatically removed the non-matching records
num_states = nrow(state_df)
state_df$state_nm = NULL
# order alphabetically by abbrev; similar to gcam-usa co2 mapping file
state_df = state_df[order(state_df$gcam_states_abr),]
state_df$out_reg_code = c(1: num_states) + num_out_reg

usa_df = merge(usa_df, state_df, by = c("state_id"), all.x = TRUE)
# this merge did not remove the non-matching records
usa_df = usa_df[usa_df$state_id <= max(state_df$state_id),]
usa_df$state_nm = NULL

max_basins = 0
mb_state_ind = NA
state_df$num_basin = NA
for(r in unique(usa_df$state_id)) {
	num_temp_basin = length(usa_df$out_basin_code[usa_df$state_id == r])
	usa_df$out_basin_code[usa_df$state_id == r] = c(1:num_temp_basin)
	state_df$num_basin[state_df$state_id == r] = num_temp_basin
	if (num_temp_basin > max_basins) {
		max_basins = num_temp_basin
		mb_state_ind = which(state_df$state_id == r)
	}

}
cat("\nMax output basins, within state", mb_state_ind, state_df$gcam_states_abr[mb_state_ind], "is", max_basins, "\n")

# set the names and order to match gcam_lu and gcam-usa co2 mapping file
usa_df = usa_df[order(usa_df$gcam_states_abr, usa_df$out_basin_code),]
names(usa_df) = c("gcam_reg_code", "gcam_basin_code", "Basin_name", "basin_abr", "land_unit", "out_lu_code", "out_basin_code", "gcam_reg_name", "out_reg_code")
usa_df = usa_df[,c("basin_abr", "gcam_reg_name", "land_unit", "out_lu_code", "gcam_reg_code", "out_reg_code", "gcam_basin_code", "Basin_name", "out_basin_code")]

# make the out lu values consistent
# but only retain the state records for further processing
gcam_state_lu = rbind(gcam_lu, usa_df)
gcam_state_lu$out_lu_code = 1:nrow(gcam_state_lu)

#try using all, cuz of USA boundaries with other countries

#gcam_state_lu = gcam_state_lu[gcam_state_lu$gcam_reg_code %in% state_df$state_id,]

# write the diagnostic output gcam region and glu codes and names

iac_gcam_state_glu_codes = gcam_state_lu[,c("out_lu_code", "land_unit","out_reg_code", "out_basin_code", "gcam_reg_name", "Basin_name", "basin_abr", "gcam_reg_code", "gcam_basin_code")]
iac_gcam_state_glu_codes = iac_gcam_state_glu_codes[order(iac_gcam_state_glu_codes$out_lu_code),]
write.table(iac_gcam_state_glu_codes, iac_state_codes_on, row.names=FALSE, quote=FALSE, sep=",")

# compile table of unique out state-regions and write the region code and region name files
# these codes may not match the actual gcam-usa region codes
# these codes are simply enumerated starting after the last region out value
gcam_state_reg_name = unique(gcam_state_lu[,c("gcam_reg_name")])
out_state_reg_list = data.frame(gcam_reg_name, stringsAsFactors=FALSE)
num_out_state_reg = nrow(out_state_reg_list)
out_state_reg_list$out_reg_code = c(1:num_out_state_reg)
write.table(out_state_reg_list$out_reg_code, state_regcode_on, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(out_state_reg_list$gcam_reg_name, state_regname_on, row.names=FALSE, col.names=FALSE, quote=FALSE)

## finish creating the usa state raster

# first filter out non-gcam states and non usa regions
in_state_rast[Which(!(in_state_rast %in% state_df$state_id), cells=TRUE)] = NA

# extract usa-only region cells from whole globe
# usa region id = 1, raster values = 1 * 10000 + glu
usa_rast = moirai_rast
usa_rast[Which(usa_rast >= 20000, cells = TRUE)] = NA

# exclude the non-usa cells from the state data
crop_state_rast = in_state_rast
crop_state_rast[Which(is.na(usa_rast), cells = TRUE)] = NA

# add the extra usa region values 
filter_rast = cover(crop_state_rast, usa_rast)

max_state_id = max(state_df$state_id)

# define a window with an odd number of pixels per side
# with equal weight 1 for each cell
s = 7
w = matrix(1,s,s)
# determine the index of the center cell
c_ind = s*(s-1)/2 + (s+1)/2

# function to replace non-state usa region value with dominant neighbor state value
fill_usa_dom <- function(x, c = c_ind) {
	# check if the focal cell (index c in nXn window) needs to be replaced with a state value
	if (x[c] > max_state_id & !is.na(x[c])) {
		# if center cell is a usa value then replace it with a neighboring state value
		# find the mode excluding NA and usa-region values
		v = x[which(!is.na(x) & x < 10000)]
		uv = unique(v)
		tab = tabulate(match(v,uv))
		return(uv[which.max(tab)])
	} else {
		return(x[c])
	}
}

# apply the function to this moving window through the raster
state_rast = focal(filter_rast, w, fun = fill_usa_dom, na.rm = FALSE, pad = TRUE)
# check the result
sp = Which(state_rast %in% state_df$state_id, cells=TRUE)
tp =  Which(!is.na(state_rast), cells=TRUE)
cat("Total state_rast pixels is", length(tp), "and state-labelled pixels is", length(sp))
cat("Target total pixels is", length(Which(!is.na(filter_rast), cells=TRUE)))
cat("Pixel # difference (132) is due to puerto rico and us virgin islands; do not include in state mapping")

# now apply the region-glu transformation to the state_rast
# this allows existing function below to work properly
# raster codes: region code * 10000 + glu
# get the basin codes
basin_rast = moirai_rast %% 10000
state_rast = state_rast * 10000 + basin_rast

# now put these state values into the whole raster so that boundaries are calculated properly
moirai_state_rast = cover(state_rast, moirai_rast)


################################################
##### first generate the gcam region map for glm

# transform moirai land unit to half-degree using dominant land unit
moirai_halfdeg_rast = aggregate(moirai_rast, fact = 6, fun = modal, na.rm = TRUE)

# extract the region values
region_halfdeg_rast = trunc(moirai_halfdeg_rast / 10000)

# convert to enumerated zones based on the gcam out xml order
out_reg_rast = region_halfdeg_rast
out_reg_rast[] = 0
for (r in 1:num_reg) {
	rcode = reg_list$gcam_reg_code[r]
	#if (rcode == taiwan_code) {
	#	# taiwan region is now china region (taiwan glus are put into China_Taiwan basin below)
	#	ocode = gcam_lu$out_reg_code[gcam_lu$gcam_reg_code == china_code][1]
	#} else {
		ocode = gcam_lu$out_reg_code[gcam_lu$gcam_reg_code == r][1]
	#}
	rinds = Which(region_halfdeg_rast == rcode, cells=TRUE)
	out_reg_rast[rinds] = ocode
}

# write to text format, then strip the header lines off of it
writeRaster(out_reg_rast, filename=reg_an, format="ascii", datatype="INT4S", overwrite=TRUE, NAflag = GLM_NODATA)
temp = read.csv(reg_an, sep=" ", skip=6, header=FALSE)
write.table(temp, reg_on, sep=" ", row.names=FALSE, col.names=FALSE)

##### second generate the gcam zone (basin/glu) map for glm

# extract the zone/glu values
glu_halfdeg_rast = moirai_halfdeg_rast %% 10000

# convert to enumerated zones based on the gcam xml order
# some moirai basins may not be in the gcam output
out_glu_rast = glu_halfdeg_rast
out_glu_rast[] = 0
for (r in 1:num_reg) {
	rcode = reg_list$gcam_reg_code[r]
	num_glu = length(moirai_codes$glu[moirai_codes$gcam_reg_code == r])
	for (g in 1:num_glu) {
		gcode = moirai_codes$glu[moirai_codes$gcam_reg_code == r][g]
		#if (rcode == taiwan_code) {
		#	ocode = gcam_lu$out_basin_code[gcam_lu$gcam_reg_code == china_code & gcam_lu$gcam_basin_code == taiwan_basin_code]
		#} else {
			ocode = gcam_lu$out_basin_code[gcam_lu$gcam_reg_code == r & gcam_lu$gcam_basin_code == gcode]
		#}
		if (length(ocode) == 0) {
			# this basin does not exist in gcam output
			ginds = Which(glu_halfdeg_rast == gcode, cells=TRUE)
			rinds = Which(region_halfdeg_rast == rcode, cells=TRUE)
			oinds = intersect(ginds, rinds)
			out_glu_rast[oinds] = 0
		} else {
			# this basin does exist in gcam output
			ginds = Which(glu_halfdeg_rast == gcode, cells=TRUE)
			rinds = Which(region_halfdeg_rast == rcode, cells=TRUE)
			oinds = intersect(ginds, rinds)
			out_glu_rast[oinds] = ocode
		}
	} # end for g loop over moirai glus
} # end for r loop over moirai regions

# write to text format, then strip the header lines off of it
writeRaster(out_glu_rast, filename=zone_an, format="ascii", datatype="INT4S", overwrite=TRUE, NAflag = GLM_NODATA)
temp = read.csv(zone_an, sep=" ", skip=6, header=FALSE)
write.table(temp, zone_on, sep=" ", row.names=FALSE, col.names=FALSE)

##### third generate the continent to region and the region to country mappings for glm
# also write a new continent file as the original doesn't match the country list length

glm_ctry = data.frame(ctrycode, ctryname)
names(glm_ctry) <- c("ctry_code", "ctry_name")
num_ctry = nrow(glm_ctry)
glm_ctry$cont_code = NA
glm_ctry$out_reg_code = NA
for (c in 1:num_ctry) {
	ccode = glm_ctry$ctry_code[c]
	cinds = Which(ctry_rast == ccode, cells=TRUE)
	cvals = ctry_rast[cinds]
	rvals = out_reg_rast[cinds]
	contvals = cont_rast[cinds]
	# use the most frequent value
	rcode = modal(rvals)
	contcode = modal(contvals)
	glm_ctry$cont_code[c] = contcode
	glm_ctry$out_reg_code[c] = rcode
	
}

# region to country file
r2c = data.frame(glm_ctry$out_reg_code, glm_ctry$ctry_code)
write.table(r2c, reg2ctry_on, row.names=FALSE, col.names=FALSE, quote=FALSE)
# new continent to ctry mapping
write.table(glm_ctry$cont_code, contcode_ctry_on, row.names=FALSE, col.names=FALSE, quote=FALSE)

# continents for out regions
# multiple answers above, so find based on region
out_reg_list$cont_code = NA
for (r in 1:num_out_reg) {
	rcode = out_reg_list$out_reg_code[r]
	rinds = Which(out_reg_rast == rcode, cells=TRUE)
	rvals = out_reg_rast[rinds]
	contvals = cont_rast[rinds]
	# use the most frequent value
	contcode = modal(contvals)
	out_reg_list$cont_code[r] = contcode
}
# new continent to region mapping
write.table(out_reg_list$cont_code, contcode_reg_on, row.names=FALSE, col.names=FALSE, quote=FALSE)

##### fourth generate the netcdf glm biomass file for gcam2glm
# copy the old file to new one so that can just replace the values
# replace lon values cuz they are off; but they are not used
system(paste0("cp ", old_bio_fn, " ", bio_on))

# now replace with the new values from bio_mat
bio = nc_open(bio_on, write=TRUE)
nbio_lon = bio$dim$lon$len
nbio_lat = bio$dim$lat$len

bio_lat = ncvar_get(bio, varid = "lat", start = c(1), count = c(nbio_lat))
bio_lon = ncvar_get(bio, varid = "lon", start = c(1), count = c(nbio_lon))
#invals = ncvar_get(bio, varid = "biomass", start = c(1,1,1), count = c(nbio_lon, nbio_lat, 1))
newlon = bio_lon + 0.25

outvals = t(bio_mat)

ncvar_put(bio, varid = "biomass", vals = outvals, start = c(1,1,1), count = c(nbio_lon, nbio_lat, 1))
ncvar_put(bio, varid = "lon", vals = newlon, start = c(1), count = c(nbio_lon))

history_out = paste0(date(), ": Alan Di Vittorio used gcam2glm_mapping.r to convert data from ", bio_fn,
	"; used the old netcdf file as a template and replaced values")
ncatt_put(bio, varid = 0, attname = "history", attval = history_out)

nc_close(bio)

##### fifth generate the gcam to glm mapping file for gcam2glm

# have to deal with the taiwan basins here as well 

# to deal with partial cells:
# get land area and cell area (km^2) of each hi-res cell on wgs84 spheroid and the cell fraction within the coarser cell
# then for each land unit calc fractions based on total weighted land unit area and total coarse cell area (which is total weighted cell area)

# half degree grid, starting at upper left corner -180,90
# cells and lon-lat limits are already aligned with the moirai data
# lat and lon are cell centers to match the initial glm data file
# so don't really have to read this in, but it is direct

# do this only if asked
if (write_gcam2glm) {

	init = nc_open(init_glm_fn)
	ninit_lon = init$dim$lon$len
	ninit_lat = init$dim$lat$len
	initlon = ncvar_get(init, varid = "lon", start = c(1), count = c(ninit_lon))
	initlat = ncvar_get(init, varid = "lat", start = c(1), count = c(ninit_lat))
	nc_close(init)

	# one degree test
	#ninit_lon = 360
	#ninit_lat = 180
	#initlon = c(-180:179) + 0.5
	#initlat = c(90:-89) - 0.5

	num_cells = ninit_lon * ninit_lat

	# run the function in parallel - this takes about 2 h 40 min for half degree
	cat("\nstarting gcam2glm mapping", date(), "\n")
	mcout = mclapply(c(0:(num_cells-1)), function(i) gcam2glm_cell(cell = i, moirai_rast, moirai_cell_area_rast, moirai_land_area_rast, gcam_lu), mc.cores = num_cores)
	cat("\nfinishing gcam2glm mapping (before sorting and writing)", date(), "\n")

	# need to merge the data frames
	gcam2glm_out = rbind.fill(mcout)

	# now sort the table by land unit and location and write it
	gcam2glm_out = gcam2glm_out[order(gcam2glm_out$out_reg_code, gcam2glm_out$out_lu_code, gcam2glm_out$center_lon, gcam2glm_out$center_lat),]
	write.table(gcam2glm_out, gcam2glm_on, row.names=FALSE, quote=FALSE, sep=",")

	cat("\nTotal output land units is", nrow(gcam_lu), "\n")

} # end if write_gcam2glm = TRUE

##### sixth generate the elm to gcam mapping file


# to deal with partial cells:
# get land area and cell (km^2) of each hi-res cell on wgs84 spheroid and the cell fraction within the coarser cell
# then for each land unit calc fractions based on total weighted land unit area and total coarse cell area (which is total weighted cell area) 

# now multiple elm resolutions are supported: 0.5 f09, and f19 (0.125 TBD)
# f09 and f19 have 0, -90 center origins
# 0.5 and0.125 have -180, -90 edge origins
# the source elm grid no longer has the cell edge values, so need to calculate them
#    need to round to address precision issues
# assume all cells are the same size (which is the case for these resolutions)

# the single cell function works for all elm resolutions as it checks for alignment with the ancillary grids before shifting/adjusting.

# match an elm surface file for 0.9424084x1.25 deg
# num lon is 288 and num lat is 192
# elm grid, edges of starting cell at -0.625, -90
# the first longitude straddles 0 (1.25 is constant size in latitude) 
# the latitude size for each pole cell is half: 0.4712042
# lon and lat increase with increasing index

# the output mapping is an index for the elm grid in the elm order

# do this only if asked
if (write_elm2gcam) {

    for(l in 1:length(clm_surf_fns)) {

		clm_surf = nc_open(clm_surf_fns[l])
	
		num_lon = clm_surf$dim$lsmlon$len
		num_lat = clm_surf$dim$lsmlat$len
		lsmlon = ncvar_get(clm_surf,varid="lsmlon",start=c(1), count=c(num_lon))
		lsmlat = ncvar_get(clm_surf,varid="lsmlat",start=c(1), count=c(num_lat))
	
		longxy = ncvar_get(clm_surf,varid="LONGXY",start=c(1,1), count=c(num_lon, num_lat))
		latixy = ncvar_get(clm_surf,varid="LATIXY",start=c(1,1), count=c(num_lon, num_lat))
	
	    nc_close(clm_surf)
	
	    # precision errors at the boundaries are checked within the single cell funciton
	
		half_lon = (longxy[2,1] - longxy[1,1]) / 2
		half_lat = (latixy[1,2] - latixy[1,1]) / 2
	
		lonE = longxy + half_lon
		lonW = longxy - half_lon
		latN = latixy + half_lat
		latS = latixy - half_lat
	
		num_cells = num_lon * num_lat

		# run the single-cell elm2gcam grid mapping function in parallel - this takes about 50 minutes
		cat("\nstarting elm2gcam mapping", l , date(), "\n")
		mcout = mclapply(c(0:(num_cells-1)), function(i) elm2gcam_cell(cell = i, moirai_rast, moirai_cell_area_rast, moirai_land_area_rast, gcam_lu), mc.cores = num_cores)
		cat("\nfinishing elm2gcam2 mapping (before sorting and writing)", l , date(), "\n")

		# need to merge the data frames
		elm2gcam_out = rbind.fill(mcout)

		# now sort the table by land unit and location and write it
		elm2gcam_out = elm2gcam_out[order(elm2gcam_out$gcam_reg_name, elm2gcam_out$xlon, elm2gcam_out$ylat, elm2gcam_out$basin_abr),]
	
		# change the column names
		names(elm2gcam_out) <- c("region_index", "GLU_index", "x", "y", "gcam_reg_name", "GLU_name", "Weight")
	
		write.table(elm2gcam_out, elm2gcam_ons[l], row.names=FALSE, quote=FALSE, sep=",")
	
		####################
		# run the single-cell function again, but with the state data
		cat("\nstarting elm2gcam usa mapping", l , date(), "\n")
		mcout = mclapply(c(0:(num_cells-1)), function(i) elm2gcam_cell(cell = i, moirai_state_rast, moirai_cell_area_rast, moirai_land_area_rast, gcam_state_lu), mc.cores = num_cores)
		cat("\nfinishing elm2gcam2 usa mapping (before sorting and writing)", l , date(), "\n")

		# need to merge the data frames
		elm2gcam_usa_out = rbind.fill(mcout)

		# now sort the table by land unit and location
		elm2gcam_usa_out = elm2gcam_usa_out[order(elm2gcam_usa_out$gcam_reg_name, elm2gcam_usa_out$xlon, elm2gcam_usa_out$ylat, elm2gcam_usa_out$basin_abr),]
	
		# change the column names
		names(elm2gcam_usa_out) <- c("region_index", "GLU_index", "x", "y", "gcam_reg_name", "GLU_name", "Weight")
	
		# check for consistency between usa region and state grid cells
		# but eliminate USA Caribbean basin from comparison because puerto rico and virgin islands are not part of the state mapping
		usa_weight = aggregate(Weight ~ x + y, data = elm2gcam_out[elm2gcam_out$gcam_reg_name == "USA" & elm2gcam_out$GLU_name != "Caribbean",], FUN=sum)
		total_state_weight = aggregate(Weight ~ x + y, data = elm2gcam_usa_out[elm2gcam_usa_out$region_index > 32,], FUN=sum)
		comp_df = merge(usa_weight, total_state_weight, by = c("x", "y"), all = TRUE)
		comp_df$diff = comp_df$Weight.x - comp_df$Weight.y
		diff_sum = sum(comp_df$diff)
		cat("sum of cell-level differences between usa region weight and total state weight is", diff_sum)
		cat("max of cell-level differences between usa region weight and total state weight is", max(comp_df$diff))
		cat("min of cell-level differences between usa region weight and total state weight is", min(comp_df$diff))
	
	    # add the state records to the region records
	    elm2gcam_usa_out = rbind(elm2gcam_out, elm2gcam_usa_out)
	
		write.table(elm2gcam_usa_out, elm2gcam_usa_ons[l], row.names=FALSE, quote=FALSE, sep=",")
	
	}

} # end if write_elm2gcam = TRUE

##### seventh, add half degree glm cell area to the initial glm file

# first read in cell_area from an original file
of = nc_open(orig_other_fn)
n_lon = of$dim$lon$len
n_lat = of$dim$lat$len
o_lon = ncvar_get(of, varid = "lon", start = c(1), count = c(n_lon))
o_lat = ncvar_get(of, varid = "lat", start = c(1), count = c(n_lat))
cell_area = ncvar_get(of, varid = "cell_area", start = c(1,1), count = c(n_lon, n_lat))
nc_close(of)

system(paste0("cp ", init_glm_fn, " ", init_glm_on))

init = nc_open(init_glm_on, write=TRUE)
ninit_lon = init$dim$lon$len
ninit_lat = init$dim$lat$len
londim = of$dim[['lon']]
latdim = of$dim[['lat']]
initlon = ncvar_get(init, varid = "lon", start = c(1), count = c(ninit_lon))
initlat = ncvar_get(init, varid = "lat", start = c(1), count = c(ninit_lat))

cell_area_var = ncvar_def("cell_area", dim = list(londim, latdim), longname = "area of grid cell", units = "m^2")
init <- ncvar_add(init, cell_area_var)
ncvar_put(init, varid = "cell_area", vals = cell_area, start = c(1,1), count = c(n_lon, n_lat))

nc_close(init)

cat("\nfinishing e3sm_gcam_land_mapping.r", date(), "\n")

} 
############ end iesmv2_land_mapping function


############# gcam2glm single-cell grid mapping function

# needs cell index and four arguments: moirai_rast, moirai_cell_area_rast, moirai_land_area_rast, gcam_lu
gcam2glm_cell <- function(cell, regglu_rast, cell_area_rast, land_area_rast, gcam_lu_df)	{
		# this assumes that cell starts at 0 to eliminate a call inside
		lt = trunc(cell / ninit_lon) + 1
		ln = cell %% ninit_lon + 1
	
		lnmin = initlon[ln] - 0.25
		lnmax = initlon[ln] + 0.25
		ltmin = initlat[lt] - 0.25
		ltmax = initlat[lt] + 0.25
		
		# one degree test
		#lnmin = initlon[ln] - 0.5
		#lnmax = initlon[ln] + 0.5
		#ltmin = initlat[lt] - 0.5
		#ltmax = initlat[lt] + 0.5
		
		# need each conrer of the larger cell, clockwise, and first and last are the same
		corners = cbind(c(lnmin, lnmax, lnmax, lnmin, lnmin), c(ltmax, ltmax, ltmin, ltmin, ltmax))
		spobj = SpatialPolygons( list(Polygons( list(Polygon(coords = corners, hole=FALSE)) , 1)) , proj4string=CRS(PROJ4_STRING))
		
		lu_cells = extract(regglu_rast, spobj, weights=TRUE, normalizeWeights=FALSE, cellnumbers=TRUE, df=TRUE, small=TRUE, exact=TRUE)
		names(lu_cells)[names(lu_cells) == names(regglu_rast)] = "gcam_lu_code"
		ca = extract(cell_area_rast, spobj, weights=TRUE, normalizeWeights=FALSE, cellnumbers=TRUE, df=TRUE, small=TRUE, exact=TRUE)
		la = extract(land_area_rast, spobj, weights=TRUE, normalizeWeights=FALSE, cellnumbers=TRUE, df=TRUE, small=TRUE, exact=TRUE)
		names(ca)[names(ca) == "layer"] = "cell_area"
		names(la)[names(la) == "layer"] = "land_area"
		lu_cells=merge(lu_cells, ca, by=c("ID", "cell"))
		lu_cells=merge(lu_cells, la, by=c("ID","cell"))
		lu_cells$weighted_land_area = lu_cells$weight.x * lu_cells$land_area
		lu_cells$weighted_cell_area = lu_cells$weight.x * lu_cells$cell_area
		# need area of entire coarse cell, but fine cells have NA where there is no land
        # so can't use this to normalize to entire grid cell!
		#total_cell_area = sum(lu_cells$weighted_cell_area, na.rm=TRUE)
		# filter out empty values and zero area
		lu_cells = lu_cells[!is.na(lu_cells$gcam_lu_code) & lu_cells$weighted_land_area > 0 & !is.na(lu_cells$weighted_land_area),]
		total_land_area = sum(lu_cells$weighted_land_area, na.rm=TRUE)
		
		#if (total_land_area > 0 & total_cell_area > 0) {
		if (total_land_area > 0) {
			
			lu_area = aggregate(weighted_land_area ~ gcam_lu_code, lu_cells, FUN=sum, na.rm=TRUE)
			lu_area$total_land_area = total_land_area
			lu_area$lu_fraction = lu_area$weighted_land_area / lu_area$total_land_area
			#lu_area$total_cell_area = total_cell_area
			#lu_area$lu_fraction = lu_area$weighted_land_area / lu_area$total_cell_area
			#lu_area$gcam_lu_code[lu_area$gcam_lu_code == taiwan_code * 10000 + taiwan_basin_code | gcam_lu_code == taiwan_code * 10000 + china_coast_basin_code] =
			#	gcam_lu_code = china_code * 10000 + taiwan_code
			weight = lu_area$lu_fraction
			gcam_reg_code = trunc(lu_area$gcam_lu_code / 10000)
			gcam_basin_code = lu_area$gcam_lu_code %% 10000
			numrec = length(weight)
			center_lon = rep(initlon[ln], numrec)
			center_lat = rep(initlat[lt], numrec)
			df = data.frame(gcam_reg_code, gcam_basin_code, center_lon, center_lat, weight)
			df = merge(df, gcam_lu_df, by=c("gcam_reg_code", "gcam_basin_code"), sort=FALSE)
			df_out = df[,c("out_reg_code", "out_lu_code", "center_lon", "center_lat", "gcam_reg_name", "basin_abr", "weight")]
			df_out$gcam_reg_name = gsub("\\s", "", df_out$gcam_reg_name)
			# if a land unit in the map does not exist in the gcam land unit list, then the land unit weights will not sum to one
			# so normalize the existing values to one
			sumw = sum(df_out$weight)
			if(sumw < 1) { df_out$weight = 1.0 / sumw * df_out$weight}
			return(df_out)
		} # end if values exist in this cell
} # end single sell gcam2glm map function


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
		
		# first check for precision errors at lat boundaries
		if (ltmin < -90) {ltmin = -90}
		if (ltmax > 90) {ltmax = 90}
		# check the lon edges for precision at -180, 180, if necessary, below
		
		# may need to shift the longitude values so that they range from -180 to 180 rather than -0.625 to 179.375
		#    deal with the boundary straddle also
		if (lnmin < 180 & lnmax > 180) {
			
			# shift if necessary
			if(lnmax > (180 + BTOL)) {
				# coarser grids
				# need two spatial objects
				# need each conrer of the larger cell, clockwise, and first and last are the same
				# lon adjust
				lnmax_tmp = lnmax
				lnmax = 180
				corners = cbind(c(lnmin, lnmax, lnmax, lnmin, lnmin), c(ltmax, ltmax, ltmin, ltmin, ltmax))
				lnmin = -180
				lnmax =  -180 + (lnmax_tmp - 180)
				corners2 = cbind(c(lnmin, lnmax, lnmax, lnmin, lnmin), c(ltmax, ltmax, ltmin, ltmin, ltmax))
				p1 = Polygons( list(Polygon(coords = corners, hole=FALSE)) , 1)
				p2 = Polygons( list(Polygon(coords = corners2, hole=FALSE)) , 2)
				spobj = SpatialPolygons( list(p1, p2) , proj4string=CRS(PROJ4_STRING))
			} else {
				# precision error on finer grid so set hard boundary for one object
				lnmax = 180
				# need each conrer of the larger cell, clockwise, and first and last are the same
				corners = cbind(c(lnmin, lnmax, lnmax, lnmin, lnmin), c(ltmax, ltmax, ltmin, ltmin, ltmax))
				spobj = SpatialPolygons( list(Polygons( list(Polygon(coords = corners, hole=FALSE)) , 1)) , proj4string=CRS(PROJ4_STRING))
			}
		} else if(lnmin < (-180 + BTOL) & lnmax > -180) {
			# this could happen only if precision error on finer grids
			lnmin = -180
			# need each conrer of the larger cell, clockwise, and first and last are the same
			corners = cbind(c(lnmin, lnmax, lnmax, lnmin, lnmin), c(ltmax, ltmax, ltmin, ltmin, ltmax))
			spobj = SpatialPolygons( list(Polygons( list(Polygon(coords = corners, hole=FALSE)) , 1)) , proj4string=CRS(PROJ4_STRING))
		} else {
			# only one spatial object needed, no precision issues but may need to shift coarser grid
			if (lnmin >= 180 & lnmax > 180) {
				lnmin = -180 + (lnmin - 180)
				lnmax = -180 + (lnmax - 180)
			}
			# need each conrer of the larger cell, clockwise, and first and last are the same
			corners = cbind(c(lnmin, lnmax, lnmax, lnmin, lnmin), c(ltmax, ltmax, ltmin, ltmin, ltmax))
			spobj = SpatialPolygons( list(Polygons( list(Polygon(coords = corners, hole=FALSE)) , 1)) , proj4string=CRS(PROJ4_STRING))
		}
		
		
		
		lu_cells = extract(regglu_rast, spobj, weights=TRUE, normalizeWeights=FALSE, cellnumbers=TRUE, df=TRUE, small=TRUE, exact=TRUE)
		names(lu_cells)[names(lu_cells) == names(regglu_rast)] = "gcam_lu_code"
		ca = extract(cell_area_rast, spobj, weights=TRUE, normalizeWeights=FALSE, cellnumbers=TRUE, df=TRUE, small=TRUE, exact=TRUE)
		la = extract(land_area_rast, spobj, weights=TRUE, normalizeWeights=FALSE, cellnumbers=TRUE, df=TRUE, small=TRUE, exact=TRUE)
		names(ca)[names(ca) == "layer"] = "cell_area"
		names(la)[names(la) == "layer"] = "land_area"
		lu_cells=merge(lu_cells, ca, by=c("ID", "cell"))
		lu_cells=merge(lu_cells, la, by=c("ID","cell"))
		lu_cells$weighted_land_area = lu_cells$weight.x * lu_cells$land_area
		lu_cells$weighted_cell_area = lu_cells$weight.x * lu_cells$cell_area
		# need area of entire coarse cell, but fine cells have NA where there is no land
		# so can't use this to normalize to entire grid cell!
		#total_cell_area = sum(lu_cells$weighted_cell_area, na.rm=TRUE)
		# filter out empty values and zero area
		lu_cells = lu_cells[!is.na(lu_cells$gcam_lu_code) & lu_cells$weighted_land_area > 0 & !is.na(lu_cells$weighted_land_area),]
		total_land_area = sum(lu_cells$weighted_land_area, na.rm=TRUE)
		
		#if (total_land_area > 0 & total_cell_area > 0) {
		if (total_land_area > 0) {
			
			lu_area = aggregate(weighted_land_area ~ gcam_lu_code, lu_cells, FUN=sum, na.rm=TRUE)
			lu_area$total_land_area = total_land_area
			lu_area$lu_fraction = lu_area$weighted_land_area / lu_area$total_land_area
			#lu_area$total_cell_area = total_cell_area
			#lu_area$lu_fraction = lu_area$weighted_land_area / lu_area$total_cell_area
			#lu_area$gcam_lu_code[lu_area$gcam_lu_code == taiwan_code * 10000 + taiwan_basin_code | gcam_lu_code == taiwan_code * 10000 + china_coast_basin_code] =
			#	gcam_lu_code = china_code * 10000 + taiwan_code
			weight = lu_area$lu_fraction
			gcam_reg_code = trunc(lu_area$gcam_lu_code / 10000)
			gcam_basin_code = lu_area$gcam_lu_code %% 10000
			numrec = length(weight)
			xlon = rep(ln, numrec)
			ylat = rep(lt, numrec)
			df = data.frame(gcam_reg_code, gcam_basin_code, xlon, ylat, weight)
			df = merge(df, gcam_lu_df, by=c("gcam_reg_code", "gcam_basin_code"), sort=FALSE)
			df_out = df[,c("out_reg_code", "out_lu_code", "xlon", "ylat", "gcam_reg_name", "basin_abr", "weight")]
			# if a land unit in the map does not exist in the gcam land unit list, then the land unit weights will not sum to one
			# so normalize the existing values to one
			sumw = sum(df_out$weight)
			if(sumw < 1) { df_out$weight = 1.0 / sumw * df_out$weight}
			return(df_out)
		} # end if values exist in this cell
} # end single sell elm2gcam map function
