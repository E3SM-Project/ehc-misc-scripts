# Author - Eva Sinha, Pacific Northwest National Lab
#
# era5_hdd_cdd_elm_grids.py
#
# Computes mean annual Heating Degree Days (HDD) and Cooling Degree Days (CDD)
# from ERA5 near-surface temperature forcing, then regrids the spatial fields
# to each of the five standard ELM grid resolutions used in E3SM-GCAM coupling.
#
# HDD/CDD method follows calc_ERA5_HDD_CDD.py:
#   - Base temperature : 18 deg C (291.15 K)
#   - HDD daily        : mean(max(0, T_ref - T_2m))  [K per day]
#   - CDD daily        : mean(max(0, T_2m - T_ref))  [K per day]
#   - Annual           : sum of daily values within each calendar year
#   - Spatial field    : mean of annual fields over YR_START:YR_END
#
# Regridding: bilinear interpolation via scipy.RegularGridInterpolator from
# ERA5 0.25 deg to each ELM grid.  For grids coarser than ERA5 (0.5, 0.9,
# 1.9 deg) this samples the interpolated ERA5 surface at each ELM cell centre;
# for finer grids (0.125 deg) it interpolates between ERA5 points.  ELM cell
# centre coordinates are read from LONGXY/LATIXY in the standard surface data
# files.  Longitudes are normalised to [-180, 180] before interpolation.
#
# WORKFLOW:
#   Step 1 - Load ERA5 hourly t2m files for YR_START:YR_END.
#            Compute per-pixel daily HDD and CDD, resample to annual sums,
#            then average over years to obtain a single (lat x lon) field.
#            Fields are loaded into memory and a RegularGridInterpolator is
#            built once for HDD and once for CDD.
#
#   Step 2 - For each ELM surface data file:
#            Read LONGXY / LATIXY to extract 1-D lon and lat cell-centre
#            arrays, build a (num_lat x num_lon) query meshgrid, interpolate
#            HDD and CDD, and write two CSV files.
#
# INPUT:
#   ERA5 t2m files : /global/cfs/cdirs/e3sm/inputdata/atm/datm7/
#                    atm_forcing.datm7.ERA.0.25d.v5.c180614/tbot/
#                    elmforc.ERA5.c2018.0.25d.t2m.YYYY-MM.nc  (hourly)
#   ELM surface    : /global/cfs/cdirs/e3sm/inputdata/iac/giac/
#                    inputs_for_mapping_generation/
#                    surfdata_<res>_*.nc  (five resolutions)
#
# OUTPUT FILES (written to ./hdd_cdd_outfiles/, year range appended to name):
#   elm0.9x1.25_HDD_YYYY_YYYY.csv,     elm0.9x1.25_CDD_YYYY_YYYY.csv
#   elm1.9x2.5_HDD_YYYY_YYYY.csv,      elm1.9x2.5_CDD_YYYY_YYYY.csv
#   elm0.5x0.5_HDD_YYYY_YYYY.csv,      elm0.5x0.5_CDD_YYYY_YYYY.csv
#   elm0.125x0.125_HDD_YYYY_YYYY.csv,  elm0.125x0.125_CDD_YYYY_YYYY.csv
#   elm0.25x0.25_HDD_YYYY_YYYY.csv,    elm0.25x0.25_CDD_YYYY_YYYY.csv
#   elm<res>_HDD_CDD_YYYY_YYYY.png     (one spatial diagnostic map per resolution)
#
# CSV columns:
#   lon_ind  - 1-based longitude index (varies fastest, matches ELM cell order)
#   lat_ind  - 1-based latitude index
#   HDD / CDD - mean annual value in °F days yr-1

from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ---------------------------------------------------------------------------
# Configuration  (mirrors calc_ERA5_HDD_CDD.py)
# ---------------------------------------------------------------------------
ERA5_TBOT_DIR = Path('/global/cfs/cdirs/e3sm/inputdata/atm/datm7/'
                     'atm_forcing.datm7.ERA.0.25d.v5.c180614/tbot/')
ELM_SURF_DIR  = Path('/global/cfs/cdirs/e3sm/inputdata/iac/giac/'
                     'inputs_for_mapping_generation')
OUT_DIR  = Path('./hdd_cdd_outfiles')

VARNAME  = 't2m'
T_REF2M  = 18.0 + 273.15   # base temperature in Kelvin (18 deg C)
YR_START = 2010
YR_END   = 2014

ELM_SURF_FNS = [
    ('surfdata_0.9x1.25_HIST_simyr2015_c201021.nc',      'base_f09_ERA5_annAvg'),
    ('surfdata_1.9x2.5_SSP5_RCP85_simyr2015_c210916.nc', 'base_f19_ERA5_annAvg'),
    ('surfdata_0.5x0.5_HIST_simyr2015_c220318.nc',       'base_r05_ERA5_annAvg'),
    ('surfdata_0.125x0.125_HIST_simyr2015_c241205.nc',   'base_r0125_ERA5_annAvg'),
    ('surfdata_0.25x0.25_simyr2015_c250312.nc',          'base_r025_ERA5_annAvg'),
]

OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Step 1 – Compute ERA5 mean-annual spatial HDD/CDD
# ---------------------------------------------------------------------------

# ----------------------------
def _normalise_lon(ds: xr.Dataset) -> xr.Dataset:
    """Shift longitudes from [0, 360] to [-180, 180] if needed."""
    if float(ds['lon'].max()) > 180:
        ds = ds.assign_coords(lon=(((ds['lon'] + 180) % 360) - 180)).sortby('lon')
    return ds

# ----------------------------
def _open_tbot(yr_start: int, yr_end: int) -> xr.Dataset:
    """Open ERA5 near-surface temperature files for the requested year range."""
    tag = 'elmforc.ERA5.c2018.0.25d.t2m.'
    years_range = set(range(yr_start, yr_end + 1))
    fpaths = sorted([
        p for p in ERA5_TBOT_DIR.glob(f'*{tag}*')
        if int(p.name.split(tag)[1][:4]) in years_range
    ])
    if not fpaths:
        raise FileNotFoundError(
            f"No ERA5 t2m files found for {yr_start}-{yr_end} in {ERA5_TBOT_DIR}")
    drop_vars = [v for v in xr.open_dataset(fpaths[0]).data_vars if v != VARNAME]
    return xr.open_mfdataset(fpaths, combine='by_coords',
                             data_vars='minimal', drop_variables=drop_vars)

# ----------------------------
def compute_hdd_cdd_spatial(yr_start: int, yr_end: int
                            ) -> tuple[xr.DataArray, xr.DataArray]:
    """Return mean annual HDD and CDD spatial fields (lat x lon, °F days yr-1).

    ERA5 t2m is in Kelvin; daily degree-days are first accumulated in K-days
    (= °C-days) and then converted to Fahrenheit-days by multiplying by 9/5,
    since a temperature difference of 1 K equals 9/5 °F.
    """
    print(f"  Loading ERA5 t2m for {yr_start}-{yr_end} ...")
    with _open_tbot(yr_start, yr_end) as ds:
        ds  = _normalise_lon(ds)
        t2m = ds[VARNAME]

        # Filter to the requested year range
        yr  = t2m.time.dt.year
        t2m = t2m.isel(time=((yr >= yr_start) & (yr <= yr_end)).values)

        print("  Computing daily HDD and CDD from hourly data ...")
        hdd_daily = (T_REF2M - t2m).clip(min=0).resample(time='D').mean()
        cdd_daily = (t2m - T_REF2M).clip(min=0).resample(time='D').mean()

        print("  Summing to annual then averaging over years ...")
        hdd_mean = hdd_daily.resample(time='1YE').sum().mean(dim='time')
        cdd_mean = cdd_daily.resample(time='1YE').sum().mean(dim='time')

        # Load into memory before closing the lazy dataset
        hdd_mean = hdd_mean.load()
        cdd_mean = cdd_mean.load()

    # Convert K-days (= °C-days) to Fahrenheit-days: ΔT_F = ΔT_K × 9/5
    hdd_mean = hdd_mean * (9.0 / 5.0)
    cdd_mean = cdd_mean * (9.0 / 5.0)

    return hdd_mean, cdd_mean
# ----------------------------

print(f"Starting ERA5 HDD/CDD on ELM grids at: {datetime.now()}\n")
print("Step 1: Computing ERA5 mean annual HDD and CDD spatial fields ...")
hdd_era5, cdd_era5 = compute_hdd_cdd_spatial(YR_START, YR_END)
print(f"  ERA5 grid: lat [{float(hdd_era5.lat.min()):.3f}, "
      f"{float(hdd_era5.lat.max()):.3f}]  "
      f"lon [{float(hdd_era5.lon.min()):.3f}, {float(hdd_era5.lon.max()):.3f}]")

# RegularGridInterpolator requires strictly ascending coordinate axes.
# ERA5 lat may be stored top-to-bottom (descending); flip if so.
if float(hdd_era5.lat[0]) > float(hdd_era5.lat[-1]):
    hdd_era5 = hdd_era5.isel(lat=slice(None, None, -1))
    cdd_era5 = cdd_era5.isel(lat=slice(None, None, -1))

era5_lats = hdd_era5.lat.values   # ascending
era5_lons = hdd_era5.lon.values   # ascending [-180, 180]

interp_hdd = RegularGridInterpolator(
    (era5_lats, era5_lons), hdd_era5.values,
    method='linear', bounds_error=False, fill_value=np.nan
)
interp_cdd = RegularGridInterpolator(
    (era5_lats, era5_lons), cdd_era5.values,
    method='linear', bounds_error=False, fill_value=np.nan
)

# ---------------------------------------------------------------------------
# Step 2 – Regrid to each ELM resolution and write CSV outputs
# ---------------------------------------------------------------------------

def plot_spatial_check(lon_1d: np.ndarray, lat_1d: np.ndarray,
                       hdd_vals: np.ndarray, cdd_vals: np.ndarray,
                       label: str, yr_start: int, yr_end: int,
                       out_path: Path) -> None:
    """Save a two-panel Cartopy map of regridded HDD and CDD for visual QA.

    Longitudes are sorted to ascending order before plotting so that
    pcolormesh renders correctly even for f09/f19 grids whose normalised
    longitudes wrap from ~180 back to ~-180.
    """
    # Sort lon to ascending order and reorder data columns to match.
    lon_sort_idx = np.argsort(lon_1d)
    lon_plot  = lon_1d[lon_sort_idx]
    hdd_plot  = hdd_vals[:, lon_sort_idx]
    cdd_plot  = cdd_vals[:, lon_sort_idx]

    # Build 2-D coordinate arrays for pcolormesh.
    lon2d, lat2d = np.meshgrid(lon_plot, lat_1d)

    fig, axes = plt.subplots(
        2, 1, figsize=(12, 8),
        subplot_kw={'projection': ccrs.PlateCarree()},
        constrained_layout=True
    )

    configs = [
        (hdd_plot, 'YlGnBu', 'HDD (°F days yr⁻¹)', 'HDD'),
        (cdd_plot, 'YlOrRd', 'CDD (°F days yr⁻¹)', 'CDD'),
    ]

    for ax, (data, cmap, cbar_label, var_title) in zip(axes, configs):
        ax.add_feature(cfeature.OCEAN, facecolor='#c9e8f0', zorder=0)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.4, edgecolor='0.35', zorder=2)
        ax.add_feature(cfeature.BORDERS,   linewidth=0.3, edgecolor='0.35', zorder=2)
        ax.set_global()

        im = ax.pcolormesh(
            lon2d, lat2d, data,
            transform=ccrs.PlateCarree(),
            cmap=cmap, vmin=0, vmax=float(np.nanmax(data)),
            shading='auto', zorder=1
        )
        plt.colorbar(im, ax=ax, orientation='horizontal',
                     pad=0.04, fraction=0.046, label=cbar_label)
        ax.set_title(f'{var_title}  {label}  ({yr_start}–{yr_end})', fontsize=11)

    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Plot:    {out_path}")


print("\nStep 2: Regridding to ELM grid resolutions ...")

for surf_fn, label in ELM_SURF_FNS:
    surf_path = ELM_SURF_DIR / surf_fn
    if not surf_path.exists():
        print(f"  Skipping {label}: surface file not found at {surf_path}")
        continue

    print(f"  Processing {label} ...")

    with xr.open_dataset(surf_path) as ds_surf:
        longxy = ds_surf['LONGXY'].values   # (lsmlat, lsmlon)
        latixy = ds_surf['LATIXY'].values
        num_lat, num_lon = longxy.shape

    # 1-D cell-centre arrays.
    # LONGXY is constant along axis 0 (rows); first row gives all lon centres.
    # LATIXY is constant along axis 1 (cols); first col gives all lat centres.
    lon_1d = longxy[0, :].copy()   # (num_lon,)
    lat_1d = latixy[:, 0].copy()   # (num_lat,)

    # Normalise ELM longitudes to [-180, 180] to match ERA5 convention.
    # f09/f19 grids use [0, 360]; 0.5/0.25/0.125 grids are already [-180, 180].
    lon_1d = np.where(lon_1d > 180.0, lon_1d - 360.0, lon_1d)
    lon_1d = np.where(lon_1d < -180.0, lon_1d + 360.0, lon_1d)

    # Build a (num_lat x num_lon) meshgrid of query points.
    # indexing='ij' keeps lat on axis 0 and lon on axis 1.
    lat_grid, lon_grid = np.meshgrid(lat_1d, lon_1d, indexing='ij')
    points = np.column_stack([lat_grid.ravel(), lon_grid.ravel()])  # (N, 2)

    hdd_vals = interp_hdd(points).reshape(num_lat, num_lon)
    cdd_vals = interp_cdd(points).reshape(num_lat, num_lon)

    # 1-based indices matching ELM cell ordering (lon varies fastest).
    lat_idx, lon_idx = np.meshgrid(
        np.arange(1, num_lat + 1),
        np.arange(1, num_lon + 1),
        indexing='ij'
    )

    yr_tag = f'{YR_START}-{YR_END}'
    for var_name, vals in [('hdd', hdd_vals), ('cdd', cdd_vals)]:
        df = pd.DataFrame({
            'lon_ind': lon_idx.ravel().astype(int),
            'lat_ind': lat_idx.ravel().astype(int),
            var_name:  vals.ravel(),
        })
        df = df.dropna(subset=[var_name])
        out_path = OUT_DIR / f'{label}_{yr_tag}_{var_name}.csv'
        df.to_csv(out_path, index=False)
        print(f"    Written: {out_path}  ({len(df)} cells)")

    plot_path = OUT_DIR / f'{label}_{yr_tag}_hdd_cdd.png'
    plot_spatial_check(lon_1d, lat_1d, hdd_vals, cdd_vals,
                       label, YR_START, YR_END, plot_path)

print(f"\nFinished at: {datetime.now()}")
