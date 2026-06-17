import netCDF4 as nc
import numpy as np
import os

__author__ = 'Eva Sinha'
__email__  = 'eva.sinha@pnnl.gov'

fpath = '/lcrc/group/e3sm/data/inputdata/lnd/clm2/paramdata/'

fname_old_25_pfts = 'clm_params_c211124.nc'
fname_new_25_pfts = 'clm_params_c211124_tuned.nc'
fname_tuned   = 'clm_params_c230517_phs_50pfts_tuned.nc'

# ---------- Read tuned parameter file ----------
ds_tuned = nc.Dataset(fpath + fname_tuned)

# ---------- Read old parameter file ----------
ds_old_25_pfts = nc.Dataset(fpath + fname_old_25_pfts)

# ---------- Create new parameter file ----------
if os.path.exists(fname_new_25_pfts):
   os.remove(fname_new_25_pfts)
    
ds = nc.Dataset(fname_new_25_pfts, 'w', format='NETCDF3_64BIT')

# ----- Copy dimensions from old parameter file -----
# https://gist.github.com/guziy/8543562 
for dname, the_dim in ds_old_25_pfts.dimensions.items():
  print (dname, len(the_dim))
  ds.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

# ----- Copy variables from old parameter file -----
for v_name, varin in ds_old_25_pfts.variables.items():
  if (v_name != 'soilpsi_off' and v_name != 'kmax'):
    print (v_name, varin.datatype, varin.dimensions) 

    # ----- Create variable  -----
    outVar = ds.createVariable(v_name, varin.datatype, varin.dimensions)

    # ----- Copy variable attributes -----
    outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    
    # ----- Copy variable values -----
    outVar[:] = varin[:]

# ----- Copy soilpsi_off and kmax from tuned parameter file -----
for v_name, varin in ds_tuned.variables.items():
  if (v_name == 'soilpsi_off'):
    outVar = ds.createVariable(v_name, varin.datatype, varin.dimensions)
    outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    outVar[:] = varin[:]

  if(v_name == 'kmax'):
    outVar = ds.createVariable(v_name, varin.datatype, varin.dimensions)
    outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    outVar[:] = varin[:, 0:25]
    
# ----- Close the output file -----
ds.close()
