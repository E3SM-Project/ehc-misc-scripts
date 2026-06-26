# Author - Eva Sinha, Pacific Northwest National Lab
#
# era5_hdd_cdd_gcam_regions.py
#
# Maps ERA5 mean-annual HDD and CDD to GCAM land units (region x GLU),
# including a second set of output files that sub-divide the USA into
# individual states following the same approach as e3sm_gcam_land_mapping.r.
#
# APPROACH:
#   The MOIRAI high-resolution land-area raster (5 arcminute, 2160x4320)
#   encodes each pixel's GCAM land unit as region_code*10000 + GLU_code.
#   ERA5 HDD/CDD fields are interpolated to that same fine grid, then a
#   land-area-weighted mean is computed for every GCAM land unit.
#   A second raster (moirai_state_rast) replaces the USA region code inside
#   the contiguous US with individual state codes, enabling state-level output.
#
# WORKFLOW:
#   Step 1 - Read GCAM context data  [done once]
#     - Parse woodharvest.xml for the ordered list of GCAM land units.
#     - Read MOIRAI_reggcam_GLU_with_spaces.csv and basin_to_country_mapping.csv.
#     - Build gcam_lu: lookup table with codes and names for every region x GLU pair.
#
#   Step 2 - Load MOIRAI raster and land-area grid  [done once]
#     - Open moirai_valid_region32_water_basin235.bsq: pixel value =
#       region_code * 10000 + GLU_code.
#     - Read hyde_land_plus.bil for land area (km^2) at each 5-arcminute pixel.
#     - Derive pixel-centre lat/lon arrays; pre-compute flat query-point array.
#
#   Step 2.5 - Build US state raster and extended land unit list  [done once]
#     - Rasterize gcamusa_state_glu_wgs84.shp (state_id + 50 offset) onto the
#       MOIRAI grid; apply a 7x7 focal-majority filter to assign a state ID to
#       every contiguous-US pixel.  Encode as state_id*10000 + GLU_code, then
#       overlay on moirai_rast to produce moirai_state_rast.
#     - Build gcam_state_lu: append one row per (state x GLU) combination to
#       gcam_lu, with state abbreviations as gcam_reg_name and out_reg_code
#       values starting after the last global region code.
#
#   Steps 3-5 repeat for each entry in PERIODS:
#
#   Step 3 - Compute ERA5 mean-annual HDD and CDD spatial fields
#     - Load ERA5 hourly t2m for yr_start:yr_end.
#     - Per pixel: daily HDD = mean(max(0, T_ref-T_2m)), CDD = mean(max(0, T_2m-T_ref))
#     - Annual sum then average over years -> one (lat x lon) field each.
#     - Build a RegularGridInterpolator for each field.
#
#   Step 4 - Interpolate ERA5 HDD/CDD to the MOIRAI 5-arcminute grid
#     - Query the interpolators at every MOIRAI pixel centre (bilinear, fill=NaN).
#
#   Step 5 - Aggregate to GCAM land units (global regions)
#     - Land-area-weighted mean per GCAM LU code; merge names; write CSV.
#
#   Step 5b - Aggregate to US state land units
#     - Same aggregation using moirai_state_rast; keep only state-level rows
#       (gcam_reg_code > 32).  Stack global rows on top; write USA states CSV.
#
# INPUT FILES:
#   /global/cfs/cdirs/e3sm/inputdata/iac/giac/inputs_for_mapping_generation/:
#     woodharvest.xml, MOIRAI_reggcam_GLU_with_spaces.csv,
#     basin_to_country_mapping.csv, moirai_valid_region32_water_basin235.bsq,
#     hyde_land_plus.bil, gcamusa_state_glu_wgs84/gcamusa_state_glu_wgs84.shp
#   /global/cfs/cdirs/e3sm/inputdata/atm/datm7/
#     atm_forcing.datm7.ERA.0.25d.v5.c180614/tbot/
#     elmforc.ERA5.c2018.0.25d.t2m.YYYY-MM.nc  (hourly)
#
# OUTPUT FILES (written to ./hdd_cdd_outfiles/, one set per period):
#   gcam_regions_HDD_CDD_<label>.csv
#     global GCAM regions: region_index, GLU_index, gcam_reg_name, GLU_name,
#                          HDD, CDD, land_area_km2
#     land_area_km2 = total MOIRAI land area for the GLU; used as weight
#     when collapsing GLUs to region level in era5_hdd_cdd_to_gcam_xml.py
#   gcam_usa_states_HDD_CDD_<label>.csv
#     global rows stacked on top of US state rows (same columns)
#   gcam_regions_HDD_moirai_<label>.png  -- ERA5 HDD on MOIRAI grid (diagnostic)
#   gcam_regions_CDD_moirai_<label>.png  -- ERA5 CDD on MOIRAI grid (diagnostic)

import os
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path

import pyproj

# Point PROJ to the active conda environment's data directory.
# Must be done before any geopandas / pyproj CRS operations.
if "CONDA_PREFIX" in os.environ:
    _proj_data = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
    if os.path.exists(_proj_data):
        os.environ["PROJ_DATA"] = _proj_data
        os.environ["PROJ_LIB"]  = _proj_data
        pyproj.datadir.set_data_dir(_proj_data)

import numpy as np
import pandas as pd
import xarray as xr
import rasterio
import geopandas as gpd
from rasterio.features import rasterize as rio_rasterize
from scipy.ndimage import generic_filter
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
INPUTDATA_DIR = Path('/global/cfs/cdirs/e3sm/inputdata')
INPUT_DIR     = INPUTDATA_DIR / 'iac/giac/inputs_for_mapping_generation'
ERA5_TBOT_DIR = (INPUTDATA_DIR / 'atm/datm7/'
                 'atm_forcing.datm7.ERA.0.25d.v5.c180614/tbot')
OUT_DIR = Path('./hdd_cdd_outfiles')

VARNAME = 't2m'
T_REF2M = 18.0 + 273.15   # base temperature in Kelvin (18 deg C)

# Each tuple: (output label, yr_start, yr_end)
PERIODS = [
    (1980, 1976, 1980),
    (1985, 1981, 1985),
    (1990, 1986, 1990),
    (1995, 1991, 1995),
    (2000, 1996, 2000),
    (2005, 2001, 2005),
    (2010, 2006, 2010),
    (2015, 2011, 2015),
]

# Full state names to abbreviations (matches R script gcam_states_names/abr)
STATE_ABBREVS = {
    "Alabama": "AL", "Alaska": "AK", "Arizona": "AZ", "Arkansas": "AR",
    "California": "CA", "Colorado": "CO", "Connecticut": "CT", "Delaware": "DE",
    "District of Columbia": "DC", "Florida": "FL", "Georgia": "GA", "Hawaii": "HI",
    "Idaho": "ID", "Illinois": "IL", "Indiana": "IN", "Iowa": "IA",
    "Kansas": "KS", "Kentucky": "KY", "Louisiana": "LA", "Maine": "ME",
    "Maryland": "MD", "Massachusetts": "MA", "Michigan": "MI", "Minnesota": "MN",
    "Mississippi": "MS", "Missouri": "MO", "Montana": "MT", "Nebraska": "NE",
    "Nevada": "NV", "New Hampshire": "NH", "New Jersey": "NJ", "New Mexico": "NM",
    "New York": "NY", "North Carolina": "NC", "North Dakota": "ND", "Ohio": "OH",
    "Oklahoma": "OK", "Oregon": "OR", "Pennsylvania": "PA", "Rhode Island": "RI",
    "South Carolina": "SC", "South Dakota": "SD", "Tennessee": "TN", "Texas": "TX",
    "Utah": "UT", "Vermont": "VT", "Virginia": "VA", "Washington": "WA",
    "West Virginia": "WV", "Wisconsin": "WI", "Wyoming": "WY",
}

OUT_DIR.mkdir(parents=True, exist_ok=True)
print(f"Starting ERA5 HDD/CDD to GCAM region mapping at: {datetime.now()}\n")

# ---------------------------------------------------------------------------
# Step 1 – Read GCAM context data  (done once)
# ---------------------------------------------------------------------------
print("Step 1: Reading GCAM context data ...")

wh_xml_fn = INPUT_DIR / 'woodharvest.xml'
tree = ET.parse(wh_xml_fn)
root = tree.getroot()
land_units = [node.text for node in root[1] if node.text is not None]

moirai_codes_fn = INPUT_DIR / 'MOIRAI_reggcam_GLU_with_spaces.csv'
moirai_codes    = pd.read_csv(moirai_codes_fn, skiprows=4)

gcam_basin2ctry_fn = INPUT_DIR / 'basin_to_country_mapping.csv'
gcam_basin2ctry    = pd.read_csv(gcam_basin2ctry_fn, skiprows=7)

reg_list = moirai_codes[['gcam_reg_code', 'gcam_reg_name']].drop_duplicates().reset_index(drop=True)

gcam_lu = pd.DataFrame({'land_unit': land_units})
gcam_lu['out_lu_code']   = range(1, len(gcam_lu) + 1)
gcam_lu['gcam_reg_name'] = None
gcam_lu['gcam_reg_code'] = None
gcam_lu['basin_abr']     = None

for _, row in reg_list.iterrows():
    mask = gcam_lu['land_unit'].str.startswith(row['gcam_reg_name'])
    gcam_lu.loc[mask, 'gcam_reg_name'] = row['gcam_reg_name']
    gcam_lu.loc[mask, 'gcam_reg_code'] = row['gcam_reg_code']
    gcam_lu.loc[mask, 'basin_abr'] = (
        gcam_lu.loc[mask, 'land_unit']
        .apply(lambda x: x[len(row['gcam_reg_name']) + 1:])
    )

out_reg_list = pd.DataFrame({'gcam_reg_name': gcam_lu['gcam_reg_name'].unique()})
out_reg_list['out_reg_code'] = range(1, len(out_reg_list) + 1)
gcam_lu = gcam_lu.merge(out_reg_list, on='gcam_reg_name', how='inner')

gcam_basin_sub = (
    gcam_basin2ctry[['GCAM_basin_ID', 'Basin_name', 'GLU_name']]
    .rename(columns={'GCAM_basin_ID': 'gcam_basin_code', 'GLU_name': 'basin_abr'})
)
gcam_lu = gcam_lu.merge(gcam_basin_sub, on='basin_abr', how='inner')
gcam_lu['out_basin_code'] = gcam_lu.groupby('out_reg_code').cumcount() + 1

num_out_reg = int(gcam_lu['out_reg_code'].max())   # = 32 global regions

print(f"  Loaded {len(gcam_lu)} GCAM land units across "
      f"{gcam_lu['gcam_reg_name'].nunique()} regions and "
      f"{gcam_lu['gcam_basin_code'].nunique()} basins.")

# ---------------------------------------------------------------------------
# Step 2 – Load MOIRAI raster and land-area grid  (done once)
# ---------------------------------------------------------------------------
print("\nStep 2: Loading MOIRAI raster and land-area grid ...")

moirai_fn = INPUT_DIR / 'moirai_valid_region32_water_basin235.bsq'

with rasterio.open(moirai_fn) as src:
    transform   = src.transform   # needed for rasterization in Step 2.5
    moirai_nrow = src.height          # 2160
    moirai_ncol = src.width           # 4320
    lon_origin  = src.bounds.left     # -180.0
    lat_origin  = src.bounds.top      #   90.0
    res_x       = src.res[0]          # 1/12 deg
    res_y       = src.res[1]          # 1/12 deg

with xr.open_dataset(moirai_fn, engine='rasterio') as ds_moirai:
    var_name    = list(ds_moirai.data_vars)[0]
    moirai_rast = ds_moirai[var_name].values[0]   # (2160, 4320), float

land_area_path   = INPUT_DIR / 'hyde_land_plus.bil'
moirai_land_area = np.fromfile(land_area_path, dtype='<f4').reshape((moirai_nrow, moirai_ncol))
moirai_land_area[moirai_land_area == -9999] = np.nan

# Pixel-centre coordinates (row 0 is northernmost; lat decreases going down).
moirai_lons = lon_origin + (np.arange(moirai_ncol) + 0.5) * res_x   # ascending
moirai_lats = lat_origin - (np.arange(moirai_nrow) + 0.5) * res_y   # descending

# Pre-compute flat query-point array for the interpolation in Step 4.
lat_grid, lon_grid = np.meshgrid(moirai_lats, moirai_lons, indexing='ij')
moirai_points = np.column_stack([lat_grid.ravel(), lon_grid.ravel()])

print(f"  MOIRAI grid: {moirai_nrow} rows x {moirai_ncol} cols  "
      f"({res_x:.4f} deg resolution)")

# ---------------------------------------------------------------------------
# Step 2.5 – Build US state raster and extended land unit list  (done once)
# ---------------------------------------------------------------------------
print("\nStep 2.5: Building US state raster and extended land unit list ...")

gcam_usa_state_fn = INPUT_DIR / 'gcamusa_state_glu_wgs84/gcamusa_state_glu_wgs84.shp'
state_gdf = gpd.read_file(gcam_usa_state_fn)
state_gdf['state_id'] = state_gdf['state_id'].astype(int) + 50   # offset to avoid region code conflicts
max_state_id = int(state_gdf['state_id'].max())

# Rasterize state polygons onto the MOIRAI grid using state_id as burn value.
in_state_arr = rio_rasterize(
    [(geom, sid) for geom, sid in zip(state_gdf.geometry, state_gdf['state_id'])],
    out_shape=(moirai_nrow, moirai_ncol),
    transform=transform,
    fill=0,
    all_touched=False,
    dtype=np.int32,
).astype(float)
in_state_arr[in_state_arr == 0] = np.nan

# 7x7 focal-majority filter: replace contiguous-US pixels that lack a state
# assignment with the most common state value among their neighbors.
def _fill_usa_dom(buf):
    center = buf[len(buf) // 2]
    if not np.isnan(center) and center > max_state_id:
        neighbors = buf[(~np.isnan(buf)) & (buf < 10000)]
        if len(neighbors) > 0:
            vals, counts = np.unique(neighbors, return_counts=True)
            return vals[np.argmax(counts)]
    return center

filter_arr    = np.where(moirai_rast >= 20000, np.nan, in_state_arr)
state_arr     = generic_filter(filter_arr, _fill_usa_dom, size=(7, 7))

basin_arr         = moirai_rast % 10000
state_arr_coded   = state_arr * 10000 + basin_arr
moirai_state_rast = np.where(np.isnan(state_arr_coded), moirai_rast, state_arr_coded)

print(f"  moirai_state_rast built. "
      f"State-coded pixels: {int(np.sum(~np.isnan(state_arr_coded)))}")

# Build gcam_state_lu: one row per (state x GLU) derived from shapefile attributes.
# Shapefile columns: state_id (already +50), state_nm (full name), glu_id, Basin_name.
state_attr = (
    state_gdf[['state_id', 'state_nm', 'glu_id']]
    .drop_duplicates()
    .rename(columns={'state_id': 'gcam_reg_code', 'glu_id': 'gcam_basin_code'})
    .copy()
)

# Map full state names to abbreviations (used as gcam_reg_name in the output).
state_attr['gcam_reg_name'] = state_attr['state_nm'].map(STATE_ABBREVS)
unmapped = state_attr['gcam_reg_name'].isna().sum()
if unmapped > 0:
    print(f"  WARNING: {unmapped} state rows could not be mapped to an abbreviation.")

# Assign out_reg_code sequentially after the 32 global regions, ordered
# alphabetically by state abbreviation (mirrors R script).
state_reg_df = (
    state_attr[['gcam_reg_code', 'gcam_reg_name']]
    .drop_duplicates()
    .sort_values('gcam_reg_name')
    .reset_index(drop=True)
)
state_reg_df['out_reg_code'] = np.arange(1, len(state_reg_df) + 1) + num_out_reg

state_attr = state_attr.merge(state_reg_df[['gcam_reg_code', 'out_reg_code']],
                               on='gcam_reg_code', how='left')

# Merge basin abbreviations from gcam_lu.
basin_lookup = gcam_lu[['gcam_basin_code', 'basin_abr']].drop_duplicates()
state_attr   = state_attr.merge(basin_lookup, on='gcam_basin_code', how='left')

# Enumerate out_basin_code within each state ordered by gcam_basin_code.
state_attr = state_attr.sort_values(['out_reg_code', 'gcam_basin_code'])
state_attr['out_basin_code'] = (
    state_attr.groupby('out_reg_code').cumcount() + 1
)

# Drop temporary column; select columns that match gcam_lu layout.
state_attr = state_attr.drop(columns='state_nm')

# Stack global + state rows and re-enumerate out_lu_code across the combined list.
_cols = ['gcam_reg_code', 'gcam_basin_code', 'gcam_reg_name',
         'basin_abr', 'out_reg_code', 'out_basin_code']
gcam_state_lu = pd.concat(
    [gcam_lu[_cols], state_attr[_cols]],
    ignore_index=True
)
gcam_state_lu['out_lu_code'] = np.arange(1, len(gcam_state_lu) + 1)

print(f"  gcam_state_lu: {len(gcam_lu)} global LUs + "
      f"{len(state_attr)} state LUs = {len(gcam_state_lu)} total")

# Subset used for the state-only aggregation merge (faster lookup).
state_lu_only = gcam_state_lu[gcam_state_lu['gcam_reg_code'] > 32].copy()

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def _normalise_lon(ds: xr.Dataset) -> xr.Dataset:
    """Shift longitudes from [0, 360] to [-180, 180] if needed."""
    if float(ds['lon'].max()) > 180:
        ds = ds.assign_coords(lon=(((ds['lon'] + 180) % 360) - 180)).sortby('lon')
    return ds


def _open_tbot(yr_start: int, yr_end: int) -> xr.Dataset:
    """Open ERA5 near-surface temperature files for the requested year range."""
    tag       = 'elmforc.ERA5.c2018.0.25d.t2m.'
    years_set = set(range(yr_start, yr_end + 1))
    fpaths    = sorted([
        p for p in ERA5_TBOT_DIR.glob(f'*{tag}*')
        if int(p.name.split(tag)[1][:4]) in years_set
    ])
    if not fpaths:
        raise FileNotFoundError(
            f"No ERA5 t2m files found for {yr_start}-{yr_end} in {ERA5_TBOT_DIR}")
    drop_vars = [v for v in xr.open_dataset(fpaths[0]).data_vars if v != VARNAME]
    return xr.open_mfdataset(fpaths, combine='by_coords',
                             data_vars='minimal', drop_variables=drop_vars)


def compute_hdd_cdd_spatial(yr_start: int, yr_end: int):
    """Return mean-annual HDD and CDD spatial fields (lat x lon, °F days yr-1).

    ERA5 t2m is in Kelvin; daily degree-days are first accumulated in K-days
    (= °C-days) and then converted to Fahrenheit-days by multiplying by 9/5,
    since a temperature difference of 1 K equals 9/5 °F.
    """
    print(f"  Loading ERA5 t2m for {yr_start}-{yr_end} ...")
    with _open_tbot(yr_start, yr_end) as ds:
        ds  = _normalise_lon(ds)
        t2m = ds[VARNAME]

        yr  = t2m.time.dt.year
        t2m = t2m.isel(time=((yr >= yr_start) & (yr <= yr_end)).values)

        print("  Computing daily HDD and CDD from hourly data ...")
        hdd_daily = (T_REF2M - t2m).clip(min=0).resample(time='D').mean()
        cdd_daily = (t2m - T_REF2M).clip(min=0).resample(time='D').mean()

        print("  Summing to annual then averaging over years ...")
        hdd_mean = hdd_daily.resample(time='1YE').sum().mean(dim='time').load()
        cdd_mean = cdd_daily.resample(time='1YE').sum().mean(dim='time').load()

    # Convert K-days (= °C-days) to Fahrenheit-days: ΔT_F = ΔT_K × 9/5
    hdd_mean = hdd_mean * (9.0 / 5.0)
    cdd_mean = cdd_mean * (9.0 / 5.0)

    return hdd_mean, cdd_mean


def _aggregate_to_lu(lu_codes_flat, hdd_flat, cdd_flat, lu_lookup):
    """Land-area-weighted mean HDD and CDD per GCAM land unit.

    Parameters
    ----------
    lu_codes_flat : 1-D array of MOIRAI raster values (region*10000 + GLU).
    hdd_flat, cdd_flat : 1-D arrays of interpolated HDD/CDD at MOIRAI pixels.
    lu_lookup : DataFrame with gcam_reg_code, gcam_basin_code, gcam_reg_name,
                basin_abr, out_reg_code, out_lu_code.

    Returns
    -------
    DataFrame with columns: region_index, GLU_index, gcam_reg_name, GLU_name,
                            HDD, CDD, land_area_km2
                            (sorted by gcam_reg_name, GLU_name).
    """
    df = pd.DataFrame({
        'gcam_lu_code': lu_codes_flat,
        'land_area':    moirai_land_area.ravel(),
        'hdd':          hdd_flat,
        'cdd':          cdd_flat,
    })
    df = df.dropna(subset=['gcam_lu_code', 'land_area', 'hdd', 'cdd'])
    df = df[df['land_area'] > 0]
    df['gcam_lu_code'] = df['gcam_lu_code'].astype(int)

    df['hdd_w'] = df['hdd'] * df['land_area']
    df['cdd_w'] = df['cdd'] * df['land_area']

    agg = df.groupby('gcam_lu_code').agg(
        total_land=('land_area', 'sum'),
        hdd_sum   =('hdd_w',    'sum'),
        cdd_sum   =('cdd_w',    'sum'),
    ).reset_index()

    agg['HDD'] = agg['hdd_sum'] / agg['total_land']
    agg['CDD'] = agg['cdd_sum'] / agg['total_land']
    agg['gcam_reg_code']   = agg['gcam_lu_code'] // 10000
    agg['gcam_basin_code'] = agg['gcam_lu_code'] % 10000

    result = agg.merge(
        lu_lookup[['gcam_reg_code', 'gcam_basin_code', 'gcam_reg_name',
                   'basin_abr', 'out_reg_code', 'out_lu_code']],
        on=['gcam_reg_code', 'gcam_basin_code'],
        how='inner',
    )
    result = result.rename(columns={
        'out_reg_code': 'region_index',
        'out_lu_code':  'GLU_index',
        'basin_abr':    'GLU_name',
        'total_land':   'land_area_km2',
    })
    out_cols = ['region_index', 'GLU_index', 'gcam_reg_name', 'GLU_name',
                'HDD', 'CDD', 'land_area_km2']
    return result[out_cols].sort_values(['gcam_reg_name', 'GLU_name']).reset_index(drop=True)


def _plot_moirai_field(data: np.ndarray, title: str, cmap: str,
                       cbar_label: str, out_path: Path) -> None:
    """Save a global Cartopy map of a field on the MOIRAI grid."""
    step      = 6   # downsample for plotting speed (~0.5 deg)
    lon_plot  = moirai_lons[::step]
    lat_plot  = moirai_lats[::step]
    data_plot = data[::step, ::step]

    if lat_plot[0] > lat_plot[-1]:   # pcolormesh needs ascending lat
        lat_plot  = lat_plot[::-1]
        data_plot = data_plot[::-1, :]

    lon2d, lat2d = np.meshgrid(lon_plot, lat_plot)

    fig, ax = plt.subplots(
        figsize=(14, 6),
        subplot_kw={'projection': ccrs.PlateCarree()},
        constrained_layout=True,
    )
    ax.add_feature(cfeature.OCEAN,     facecolor='#c9e8f0', zorder=0)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.4, edgecolor='0.35', zorder=2)
    ax.add_feature(cfeature.BORDERS,   linewidth=0.3, edgecolor='0.35', zorder=2)
    ax.set_global()

    im = ax.pcolormesh(
        lon2d, lat2d, data_plot,
        transform=ccrs.PlateCarree(),
        cmap=cmap, vmin=0, vmax=float(np.nanmax(data_plot)),
        shading='auto', zorder=1,
    )
    plt.colorbar(im, ax=ax, orientation='horizontal',
                 pad=0.04, fraction=0.046, label=cbar_label)
    ax.set_title(title, fontsize=11)
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Plot:    {out_path}")

# ---------------------------------------------------------------------------
# Main loop – Steps 3-5b repeated for each time period
# ---------------------------------------------------------------------------
print(f"\nProcessing {len(PERIODS)} time periods ...\n")

for label, yr_start, yr_end in PERIODS:
    print(f"{'='*60}")
    print(f"Period: {label}  ({yr_start}–{yr_end})")

    # --- Step 3: ERA5 HDD/CDD spatial fields ---
    hdd_era5, cdd_era5 = compute_hdd_cdd_spatial(yr_start, yr_end)
    print(f"  ERA5 grid: lat [{float(hdd_era5.lat.min()):.3f}, "
          f"{float(hdd_era5.lat.max()):.3f}]  "
          f"lon [{float(hdd_era5.lon.min()):.3f}, {float(hdd_era5.lon.max()):.3f}]")

    if float(hdd_era5.lat[0]) > float(hdd_era5.lat[-1]):   # ensure ascending lat
        hdd_era5 = hdd_era5.isel(lat=slice(None, None, -1))
        cdd_era5 = cdd_era5.isel(lat=slice(None, None, -1))

    interp_hdd = RegularGridInterpolator(
        (hdd_era5.lat.values, hdd_era5.lon.values), hdd_era5.values,
        method='linear', bounds_error=False, fill_value=np.nan,
    )
    interp_cdd = RegularGridInterpolator(
        (cdd_era5.lat.values, cdd_era5.lon.values), cdd_era5.values,
        method='linear', bounds_error=False, fill_value=np.nan,
    )

    # --- Step 4: Interpolate to MOIRAI grid ---
    print("  Step 4: Interpolating to MOIRAI 5-arcminute grid ...")
    hdd_moirai = interp_hdd(moirai_points).reshape(moirai_nrow, moirai_ncol)
    cdd_moirai = interp_cdd(moirai_points).reshape(moirai_nrow, moirai_ncol)
    print(f"  HDD range: [{np.nanmin(hdd_moirai):.1f}, {np.nanmax(hdd_moirai):.1f}] °F days yr-1")
    print(f"  CDD range: [{np.nanmin(cdd_moirai):.1f}, {np.nanmax(cdd_moirai):.1f}] °F days yr-1")

    hdd_flat = hdd_moirai.ravel()
    cdd_flat = cdd_moirai.ravel()

    # --- Step 5: Global aggregation ---
    print("  Step 5: Aggregating to global GCAM land units ...")
    result = _aggregate_to_lu(moirai_rast.ravel(), hdd_flat, cdd_flat, gcam_lu)

    csv_path = OUT_DIR / f'gcam_regions_HDD_CDD_{label}.csv'
    result.to_csv(csv_path, index=False)
    print(f"  Written {len(result)} global land units to: {csv_path}")

    # --- Step 5b: US state aggregation ---
    print("  Step 5b: Aggregating to US state land units ...")

    state_result = _aggregate_to_lu(
        moirai_state_rast.ravel(), hdd_flat, cdd_flat, state_lu_only
    )

    # Stack global rows on top of state rows (mirrors R script rbind order).
    usa_result = pd.concat([result, state_result], ignore_index=True)

    usa_csv_path = OUT_DIR / f'gcam_usa_states_HDD_CDD_{label}.csv'
    usa_result.to_csv(usa_csv_path, index=False)
    print(f"  Written {len(result)} global + {len(state_result)} state rows "
          f"to: {usa_csv_path}")

    # --- Diagnostic maps ---
    _plot_moirai_field(
        hdd_moirai,
        title=f'ERA5 HDD interpolated to MOIRAI grid  ({yr_start}–{yr_end})',
        cmap='YlGnBu', cbar_label='HDD (°F days yr⁻¹)',
        out_path=OUT_DIR / f'gcam_regions_HDD_moirai_{label}.png',
    )
    _plot_moirai_field(
        cdd_moirai,
        title=f'ERA5 CDD interpolated to MOIRAI grid  ({yr_start}–{yr_end})',
        cmap='YlOrRd', cbar_label='CDD (°F days yr⁻¹)',
        out_path=OUT_DIR / f'gcam_regions_CDD_moirai_{label}.png',
    )

print(f"\nFinished at: {datetime.now()}")
