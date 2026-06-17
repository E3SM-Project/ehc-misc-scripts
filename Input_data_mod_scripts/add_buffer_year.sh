# Copy original, create unlimited dimension
drc_in=/lcrc/group/e3sm/public_html/inputdata/atm/cam/ggas
fl_in=ne30pg2_CSEM_historical_ocean_flux_1849-2014_c20240225.nc
fl_out=ne30pg2_CSEM_historical_ocean_flux_1849-2015_c20260128.nc

/bin/cp ${drc_in}/${fl_in} ~/tmp.nc
ncks -O --mk_rec_dim time ~/tmp.nc ~/tmp.nc

# Create 2015 from 2014
ncrcat -O -d time,-12,-1 ~/tmp.nc ~/tmp_2015.nc
ncatted -a units,time,o,c,"days since 0002-01-01 00:00:00" ~/tmp_2015.nc
ncap2 -O -s 'date+=10000' ~/tmp_2015.nc ~/tmp_2015.nc

# Concatenate 2015
ncrcat -O  ~/tmp.nc ~/tmp_2015.nc ~/${fl_out}
