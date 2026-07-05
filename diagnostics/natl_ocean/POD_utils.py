# Import Analysis Tools
import os
import json
from pathlib import Path
import numpy as np
import xarray as xr
import xesmf as xe
from scipy import stats
import gsw_xarray as gsw
from numba import guvectorize
import cftime
import xwmt
import pandas as pd

# AMOC analysis: xhistogram for density-coordinate binning, cmip_basins for
# Atlantic / Indo-Pacific / Global basin masks.
from xhistogram.xarray import histogram
import cmip_basins

# Import Plotting Tools
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
from matplotlib.colors import BoundaryNorm

# Import Colors
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap

# Warnings are hidden with the below code. Comment out if you want warnings
import warnings
warnings.filterwarnings('ignore')

# Import additional utilities
import spna_masks


##### DATA LOADING (native / non-CMORized) #####
def load_native_timeslice(catalog_json, varnames, startdate, enddate):
    """Load native POP (CESM timeslice) files and harmonize to CMIP variables/coords.

    CESM timeslice output is not CMORized: variables carry native POP names (TEMP,
    SALT, UVEL, ...) on native coords (z_t in cm, TLONG/TLAT, ULONG/ULAT), with native
    units, and there is no tos. ESNB opens by CMIP variable_id and requires that name
    to exist in-file, so we bypass it here and harmonize in-memory. Shared by both the
    ESNB and intake-esm notebooks so the two stay consistent.

    The native<->CMIP mapping is read from the catalog's `native_variable_name` column
    (data-driven). The resulting `mapping` dict has the same {cmip: native} shape as
    e.g. `mapping = {"zg": "Z500"}`.

    Harmonization steps (CESM/POP only -- NOT applied to CMORized models):
      * rename native variable -> CMIP id and native coords -> lon/lat/lev,
      * normalize units (so g/kg, hfds W m-2, wfo kg m-2 s-1, uo/vo cm/s -> m/s),
      * recenter POP end-of-interval monthly time stamps to the interval midpoint,
      * NCAR/POP2 surface-flux recipe (mirrors sub2sub_HR/wmt.py + the OSNAP FOSI
        recipe): hfds = SHF + QFLUX and wfo = SFWF - QFLUX/(latent_heat_fusion/1e4),
        where the catalog maps SHF->hfds, QFLUX->hfsifrazil, SFWF->wfo,
      * derive tos from the surface level of thetao (native output has no SST field).

    Returns a list of single-variable xr.Datasets (one per CMIP variable) that the
    per-variable dispatch consumes. preprocess_coords / check_depth_units finish
    coordinate harmonization downstream (nlat/nlon->y/x, lev cm->m).
    """
    # POP2 latent heat of fusion (read from CESM FOSI output: 3.337e9 erg/g).
    # /1e4 converts erg/g -> J/kg, matching the OSNAP NCAR frazil recipe.
    LATENT_HEAT_FUSION_ERG_PER_G = 3.337e9

    # Resolve the catalog CSV path from the catalog JSON metadata.
    with open(catalog_json) as f:
        meta = json.load(f)
    csv_path = meta['catalog_file']
    if csv_path.startswith('file://'):
        csv_path = csv_path[len('file://'):]
    if not os.path.isabs(csv_path):
        csv_path = str(Path(catalog_json).parent / csv_path)
    df = pd.read_csv(csv_path)

    # Build the native<->CMIP mapping from the catalog column:
    #   mapping = {cmip_variable_id: native_variable_name}, e.g. {"thetao": "TEMP"}
    mapping = dict(zip(df['variable_id'], df['native_variable_name']))
    print('mapping (CMIP -> native):', mapping)

    _coord_rename = {'TLONG': 'lon', 'TLAT': 'lat',
                     'ULONG': 'lon', 'ULAT': 'lat', 'z_t': 'lev'}
    _unit_attr_fix = {'gram/kilogram': 'g/kg',        # SALT (== psu numerically)
                      'Watts/meter^2': 'W m-2',       # QFLUX
                      'watt/m^2': 'W m-2',            # SHF
                      'kg/m^2/s': 'kg m-2 s-1'}        # SFWF

    # Load every variable in the catalog (includes hfsifrazil, needed by the recipe).
    by_var = {}
    for var in mapping:
        rows = df[df['variable_id'] == var]
        if rows.empty:
            continue
        native = mapping[var]
        ds = xr.open_dataset(rows.iloc[0]['path'], chunks={'time': 1})  # lazy/dask
        ds = ds.rename({native: var})                       # native -> CMIP variable name
        ds = ds.rename({k: v for k, v in _coord_rename.items() if k in ds.variables})
        # Normalize native POP units to CMIP conventions. Velocities convert
        # cm/s -> m/s (factor 0.01); the rest are attribute-only fixes (no value change).
        if var in ('uo', 'vo') and ds[var].attrs.get('units') in ('centimeter/s', 'cm/s'):
            ds[var] = (ds[var] * 0.01).assign_attrs({**ds[var].attrs, 'units': 'm s-1'})
        _u = ds[var].attrs.get('units')
        if _u in _unit_attr_fix:
            ds[var].attrs['units'] = _unit_attr_fix[_u]
        # POP stamps monthly time at the END of the averaging interval; recenter to the
        # interval midpoint (from time_bound) so date-range selection is intuitive.
        if 'time_bound' in ds.variables:
            ds = ds.assign_coords(time=ds['time_bound'].mean(ds['time_bound'].dims[-1]))
        ds = ds.sel(time=slice(startdate, enddate))
        by_var[var] = ds[[var]]

    # NCAR/POP2 surface-flux recipe (CESM-only; mirrors the OSNAP FOSI NCAR recipe).
    # hfsifrazil (QFLUX) is the frazil ice-formation heat flux. It is added back to the
    # net surface heat flux (SHF) to form the CMIP hfds, and converted to an equivalent
    # freshwater flux (QFLUX / (latent_heat_fusion/1e4), in kg m-2 s-1) that is removed
    # from the freshwater flux (SFWF) to form wfo.
    if 'hfsifrazil' in by_var:
        frazil = by_var['hfsifrazil']['hfsifrazil']
        if 'hfds' in by_var:
            hfds = by_var['hfds']['hfds']
            by_var['hfds']['hfds'] = (hfds + frazil).assign_attrs(hfds.attrs)
        if 'wfo' in by_var:
            latfus = LATENT_HEAT_FUSION_ERG_PER_G / 1.0e4   # erg/g -> J/kg
            wfo = by_var['wfo']['wfo']
            by_var['wfo']['wfo'] = (wfo - frazil / latfus).assign_attrs(wfo.attrs)
        del by_var['hfsifrazil']   # consumed by the recipe; not a POD input

    datasets = list(by_var.values())

    # Derive tos (SST) from the surface level of thetao -- native output has no tos.
    thetao_ds = next((d for d in datasets if 'thetao' in d.data_vars), None)
    if thetao_ds is not None:
        tos = thetao_ds['thetao'].isel(lev=0, drop=True).rename('tos').to_dataset()
        datasets.append(tos)

    return datasets


##### PART 1 PROCESSING  #####
def preprocess_coords(ds):
    rename_coords_dict = {
        'd2': 'bnds',
        'axis_nbounds': 'bnds',
        'olevel': 'lev',
        'olevel_bounds': 'lev_bnds',
        'lev_bounds': 'lev_bnds',
        'time_bounds': 'time_bnds',
        'bounds_lon': 'lon_bnds',
        'bounds_lat': 'lat_bnds',
        'bounds_nav_lon': 'lon_bnds',
        'bounds_nav_lat': 'lat_bnds',
        'nav_lon': 'lon',
        'nav_lat': 'lat',
        'latitude': 'lat',
        'longitude': 'lon',
        'longitude_bnds':'lon_bnds',
        'latitude_bnds':'lat_bnds',
        'i': 'x',
        'j': 'y',
        'nlat':'y',
        'nlon':'x',
        'z_l': 'lev',
        'z_i': 'lev_bnds',
        'xh': 'lon',
        'yh': 'lat',
    }
    
    # rename when possible
    for k, v in rename_coords_dict.items():
        try:
            ds = ds.rename({k: v})
        except:
            pass
            
    ## May not be needed in general (just for CESM?)
    ds = check_depth_units(ds)

    return ds

def compute_sigma0(da_t,da_s):
    """
    Compute potential density anomaly (sigma0) from conservative temperature and salinity.

    Parameters:
    da_t : xarray.DataArray
        Conservative temperature in degrees Celsius ('degC').
    da_s : xarray.DataArray
        Absolute salinity in g/kg. Accepts 'psu' or '0.001' as equivalent.

    Returns:
    xarray.DataArray
        Sigma0 density anomaly (kg/m³ - 1000).

    Raises:
    Exception
        If input units are not recognized.
    """
    if (da_t.attrs['units']!='degC'):
            raise Exception("Check units of thetao!")
    if (da_s.attrs['units']!='g/kg'):
        if (da_s.attrs['units']=='psu') or (da_s.attrs['units']=='0.001'):
            da_s = da_s.assign_attrs({'units':'g/kg'})
        else:
            raise Exception("Check units of so!")
    da_sig0 = gsw.density.sigma0(SA=da_s,CT=da_t)
    return da_sig0


def compute_mld(da_sig0, dsig=0.03, rdep=10.0):
    """
    Compute mixed layer depth (MLD) based on a density threshold criterion.

    Parameters:
    da_sig0 : xarray.DataArray
        Potential density (sigma0) with units 'kg/m^3'.
    dsig : float, optional
        Density increase threshold from reference depth (default: 0.03 kg/m³).
    rdep : float, optional
        Reference depth in meters (default: 10.0 m).

    Returns:
    xarray.DataArray
        Mixed layer depth in meters.

    Raises:
    Exception
        If input units are not 'kg/m^3'.
    """
    if (da_sig0.attrs['units']!='kg/m^3'):
        raise Exception("Check units of sigma0!")
    zdim = 'lev'
    da_sig0_ref = da_sig0.interp(lev=rdep)
    z=da_sig0[zdim]
    k_ref = int(np.min(np.argwhere(z.values >= rdep)))
    da_mld = xr.apply_ufunc(mld_gufunc, da_sig0, z, k_ref, da_sig0_ref, dsig, 
                            input_core_dims=[[zdim], [zdim], [], [], []], dask='parallelized')
    return da_mld


@guvectorize(['void(float64[:], float64[:], intp, float64,  float64,float64[:])'], '(n),(n),(),(),()->()', nopython=True)
def mld_gufunc(phi,  z, k_ref, phi_ref,  phi_step, res):
    #  function from Guillaume Serazin via Anne-Marie Treguier
    #  phi is an xarray of density sigma0
    #  z is the depth
    #  k_ref is the reference level (first index below the reference depth)
    #  phi_ref is the density at the reference level 
    #  phi_step=0.03 (default) is the density threshold to define the mixed layer depth
    #  res is the result = mixed layer depth    
    phi_mlb = phi_ref + phi_step
    if (np.isfinite(phi[k_ref])) :
        jkmax = np.max(np.argwhere(np.isfinite(phi)))
        #  find the indices of points below kref where tabprof and reference 
        #  differ by more than zdel (use absolute value)
        indices=np.argwhere((phi[k_ref:] >= phi_mlb))
        if (indices.size == 0) :
            #  no values differ by more than zdel: mldepth is depth of the last level
            res[0]=z[jkmax]
        else :
            k_mlb = np.min(indices)+k_ref   
            # Make a linear interpolation to find the approximate MLD
            # the logic is to find a point on a 
            delta_z = z[k_mlb - 1] - z[k_mlb]
            alpha = delta_z/(phi[k_mlb - 1] - phi[k_mlb])
            beta = z[k_mlb] - alpha * phi[k_mlb]
            res[0] = alpha * phi_mlb + beta
    else :
        res[0] = np.nan


def get_dlev(lev, lev_bnds, depth_limit=10000):
    """
    Compute vertical layer thickness from depth bounds, with optional depth clipping.

    Parameters:
    lev_bnds : xarray.DataArray
        Depth bounds with shape (..., 2), where the last or first dimension corresponds to layer bounds.
    depth_limit : float, optional
        Maximum depth for clipping bounds (default: 10000 m).

    Returns:
    xarray.DataArray
        Layer thickness (dlev) in meters.
    """
    lev_bnds_lim = lev_bnds.where(lev_bnds < depth_limit, depth_limit)
    if (np.size(lev_bnds_lim.dims)==1):
        dlev = xr.DataArray(lev_bnds_lim[1:].values - lev_bnds_lim[0:-1].values,coords={'lev':lev})
    elif (np.size(lev_bnds_lim.dims)==2):
        dlev = lev_bnds_lim.isel(bnds=1) - lev_bnds_lim.isel(bnds=0)
    elif (np.size(lev_bnds_lim.dims)==3):
        dlev = xr.DataArray(lev_bnds_lim[1:].values - lev_bnds_lim[0:-1].values) #,coords={'lev':lev})
    else:
        raise ValueError('ERROR: could not handle lev_bnds')
    return dlev

def check_depth_units(ds):
    # Convert lev if needed
    if 'lev' in ds.coords:
        lev = ds['lev']
        if lev.max() > 8000:
            lev_converted = lev / 100
            lev_converted.attrs.update(lev.attrs)
            lev_converted.attrs['units'] = 'm'
            ds = ds.assign_coords(lev=lev_converted)

    # Convert lev_bnds if needed
    if 'lev_bnds' in ds.variables:
        lev_bnds = ds['lev_bnds']
        if lev_bnds.max() > 8000:
            lev_bnds_converted = lev_bnds / 100
            lev_bnds_converted.attrs.update(lev_bnds.attrs)
            lev_bnds_converted.attrs['units'] = 'm'
            ds['lev_bnds'] = lev_bnds_converted

    return ds

def compute_zavg(ds, var, dz, depth=200):
    """
    Compute thickness-weighted mean over depth for field.
    """
    # Ensure lev is available
    lev = ds['lev']
    # Find valid levels shallower than depth
    valid_mask = lev <= depth
    valid_levs = lev.where(valid_mask, drop=True)
    if valid_levs.size == 0:
        print(f"Warning: No valid levels found for depth={depth}m in variable '{var}'")
        return ds

    # Slice the variable and weights to these levels. dz is derived from the fx
    # grid (volcello/areacello) while the tracer field comes from the ESNB-loaded
    # monthly grid; their lev coords can differ by float-rounding even after the
    # cm->m conversion, which breaks label-based dz.sel. Select dz by nearest
    # level, then force its lev coord to the data's so the weighted mean aligns.
    data = ds[var].sel(lev=valid_levs)
    dz_sel = dz.sel(lev=valid_levs, method='nearest')
    dz_sel = dz_sel.assign_coords(lev=data['lev'])
    dz_sel = dz_sel.rename({'nlat': 'y', 'nlon': 'x'})

    # Do the weighted mean
    newvar = f"{var}_zavg"

    # Add to dataset and annotate
    ds[newvar] = data.weighted(dz_sel).mean('lev', keep_attrs=True).astype('float32')
    ds[newvar].attrs['zavg'] = f'0-{depth}m'

    return ds


def regrid(ds, var=None, names=None, target_grid=None, dlon=1, dlat=1, method='bilinear'):
    """
    Regrid a dataset or variable to a regular lat-lon grid using xESMF.

    Parameters:
    ds : xarray.Dataset
        Input dataset containing 'lat' and 'lon' coordinates.
    var : str, optional
        Name of the variable to regrid. If None, regrids the entire dataset.
    names : unused
        Placeholder parameter (not currently used).
    target_grid : dict or xarray.Dataset, optional
        Target grid definition. If None, creates a global grid with given resolution.
    dlon : float, optional
        Longitudinal resolution for default grid (default: 1°).
    dlat : float, optional
        Latitudinal resolution for default grid (default: 1°).
    method : str, optional
        Regridding method (e.g., 'bilinear', 'nearest_s2d', 'conservative').

    Returns:
    xarray.DataArray or xarray.Dataset
        Regridded variable or dataset.

    Notes:
    - Filters invalid lat/lon values before regridding.
    - Uses periodic boundary conditions and ignores degenerate grids.
    """
    if target_grid is None:
        target_grid = xe.util.grid_global(dlon, dlat, cf=True, lon1=360)
    
    # ds = utils.regrid(self[k], var=vn, target=target_grid, method=method)
    dsgrid = ds[['lat', 'lon']].copy()

    ## some grids (GFDL,NCAR) have large fill values that should actually be NaN
    lat = dsgrid['lat']
    lon = dsgrid['lon']
    dsgrid['lat'] = lat.where(lat>-90).where(lat<90)
    dsgrid['lon'] = lon.where(lon>-360).where(lon<360)

    regridder = xe.Regridder(dsgrid, target_grid, method=method, periodic=True, ignore_degenerate=True)
    if var is None:
        rgd = regridder(ds, keep_attrs=True, skipna=True, na_thres=0.6)
    else:
        rgd = regridder(ds[var], keep_attrs=True, skipna=True, na_thres=0.6)

    return rgd


def forcing_cycles(expid,nt):
    ''' Function to determine the number and year range of forcing
    cycles included in this OMIP data array.'''
    maxcycles = 10
    if (expid=="omip1" or expid=="omip1-spunup"):
        # For OMIP1, we expect a 1948-2009 (62-year) forcing cycle = 744 mon
        # but some groups submitted a 1948-2007 (60-year) cycle = 720 mon
        # and MIROC submitted a 1958-2009 (52 year) cycle = 624 mon
        if (nt % 744 == 0):
            nyear = 62
            yearrange = (1948,2009)
        elif (nt % 720 == 0):
            nyear = 60
            yearrange = (1948,2007)
        elif (nt % 624 == 0):
            nyear = 52
            yearrange = (1958,2009)
        else:
            raise ValueError('ERROR: could not determine OMIP1 forcing cycle')
    elif (expid=="omip2" or expid=="omip2-spunup"):
        # For OMIP2, we expect a 1958-2018 (61-year) forcing cycle = 732 mon
        if (nt % 732 == 0):
            nyear = 61
            yearrange = (1958,2018)
        elif (nt % 624 == 0):
            nyear = 52
            yearrange = (1958,2009)
        else:
            raise ValueError('ERROR: could not determine OMIP2 forcing cycle')
    else:
        raise ValueError('ERROR: experiment_id not recognized')
    ntmon = nyear*12*np.arange(1,maxcycles+1,1)
    findcyc = np.where(nt==ntmon)
    ncyc = findcyc[0][0] + 1
    print('  found '+str(ncyc)+' forcing cyles spanning '+str(yearrange))
    return ncyc,yearrange


##### PART 1 PLOTTING #####
def blue2red_cmap(n):
    """ combine two existing color maps to create a diverging color map with white in the middle
    n = the number of contour intervals
    """
    if (int(n/2) == n/2):
        # even number of contours
        nwhite=1
        nneg=n/2
        npos=n/2
    else:
        nwhite=2
        nneg = (n-1)/2
        npos = (n-1)/2
    colors1 = plt.cm.Blues_r(np.linspace(0,1, int(nneg)))
    colors2 = plt.cm.YlOrRd(np.linspace(0,1, int(npos)))
    colorsw = np.ones((nwhite,4))
    colors = np.vstack((colors1, colorsw, colors2))
    cmap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    return cmap

def get_varname(var):
    vardict = {'thetao_zavg':'T200','so_zavg':'S200','sigma0_zavg':r'${\sigma_0}$200','mld':'MLD'}
    return vardict[var]

def SpatialBias_panel(da, stats, focus_region, focus_model, ax):
    var = da.name
    var_name = get_varname(var)
    da_model = da.sel(model=focus_model)  #Create DataArray of the target model alone
    rmse_model = stats[var+'_rmse'].sel(model=focus_model).values
    bias_model = stats[var+'_bias'].sel(model=focus_model).values

    # Fontsizes
    font_title = 18
    font_label = 16
    font_tick = 16
    font_stat = 14

    # Settings
    proj = ccrs.PlateCarree()
    facecolor="grey"
        
    # set up contour levels and color map
    if var=='so_zavg':
        ci = 0.25; cmin = -2.5; cmax = 2.5
    elif var=='thetao_zavg':
        ci = 0.5; cmin = -8; cmax = 8
    elif var=='sigma0_zavg':
        ci = 0.1; cmin = -2; cmax = 2
    else:
        ci = 50; cmin = -500; cmax = 500
    nlevs = (cmax-cmin)/ci + 1
    levs = np.arange(cmin, cmax+ci, ci)
    cmap = cmap = blue2red_cmap(nlevs)
    cmap.set_over('magenta')
    cmap.set_under('cyan')
    norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)

    # Setup
    ax.set_aspect('auto')
    ax.set_facecolor(facecolor)
    ax.add_feature(cfeature.COASTLINE)

    # Plot
    cntr1 = ax.pcolormesh(da_model['lon'], da_model['lat'], da_model, shading='nearest', cmap=cmap, norm=norm, rasterized=True, transform=proj)
    cbar = plt.colorbar(cntr1, ax=ax, orientation='vertical', spacing='proportional', pad=.01)
    cbar.set_label(label=r'{}'.format(da.units), size=font_label, rotation=270, labelpad=5)
    cbar.ax.tick_params(labelsize=font_tick)

    # Show focus region
    r = focus_region
    ax.add_patch(mpatches.Rectangle(xy=[r[0], r[2]], width=abs(r[0]-r[1]), height=abs(r[3]-r[2]),
                                    facecolor=None, edgecolor='lime', fill=False, linewidth=2,
                                    transform=ccrs.PlateCarree()))

    gl = ax.gridlines(draw_labels=True, dms=True,  linewidth=0.5, color='k', alpha=0.8)
    gl.xlabel_style = {'fontsize': font_tick}
    gl.ylabel_style = {'fontsize': font_tick}
    gl.top_labels = False
    gl.right_labels = False
    ax.set_title('{}:\nClimatological {} Bias ({})'.format(focus_model,var_name,da.time_avg), fontsize=font_title)
    ax.set_title(r'rmse={:.1f}'.format(rmse_model) + '\n' + r'bias={:.1f}'.format(bias_model), fontsize=font_stat, loc='right')

def SpatialRank_panel(da, focus_region, focus_model, ax):
    var = da.name
    var_name = get_varname(var)
    da_model = da.sel(model=focus_model)  # Create DataArray of the target model alone

    # Fontsizes
    font_title = 18
    font_label = 16
    font_tick = 16

    # Settings
    proj = ccrs.PlateCarree()
    facecolor="grey"
    levs = np.arange(0, 110, 10)  # Colorbar Levels
    nlevs = levs.size
    cmap = cmap = blue2red_cmap(nlevs)
    norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)

    # Setup
    ax.set_aspect('auto')
    ax.set_facecolor(facecolor)
    ax.add_feature(cfeature.COASTLINE)

    # Process Rank Data
    rank_da = xr.apply_ufunc(stats.percentileofscore, # I think default is kind='rank' & nan_policy='propogate'
                             abs(da.where(da.model!=focus_model, drop=True)),
                             abs(da_model),
                             'rank',
                             'omit',
                             input_core_dims=[['model'], [], [], []],
                             vectorize=True,
                            )

    # Plot
    cntr1 = ax.pcolormesh(rank_da['lon'], rank_da['lat'], rank_da, shading='nearest', cmap=cmap, norm=norm, rasterized=True, transform=proj)
    cbar = plt.colorbar(cntr1, ax=ax, orientation='vertical', spacing='proportional', pad=0.01)
    cbar.set_label(label='%', size=font_label, rotation=270, labelpad=5)
    cbar.ax.tick_params(labelsize=font_tick)

    # Show focus region
    r = focus_region
    ax.add_patch(mpatches.Rectangle(xy=[r[0], r[2]], width=abs(r[0]-r[1]), height=abs(r[3]-r[2]),
                                    facecolor=None, edgecolor='lime', fill=False, linewidth=2,
                                    transform=ccrs.PlateCarree()))

    # Add lat/lon grid
    gl = ax.gridlines(draw_labels=True, dms=True,  linewidth=0.5, color='k', alpha=0.8)
    gl.xlabel_style = {'fontsize': font_tick}
    gl.ylabel_style = {'fontsize': font_tick}
    gl.top_labels = False
    gl.right_labels = False
    ax.set_title('{}:\nClimatological {} Bias Rank ({})'.format(focus_model,var_name,da.time_avg), fontsize=font_title)

def Scatter_panel(da1, da2, focus_model, ax):
    # Fontsizes
    font_title = 14
    font_label = 12
    font_tick = 12
    font_stat = 12

    # Strip out focus model
    da1_foc = da1.sel(model=focus_model)
    da2_foc = da2.sel(model=focus_model)
    da1 = da1.drop_sel(model=focus_model)
    da2 = da2.drop_sel(model=focus_model)

    mymarkers=[".","o","v","^","<",">","s","p","P","h","H","+","x","X","D","d","1","2"]
    mycolors = ['r','b','g','m']

    for i, key in enumerate(da1.model):  # Add each model to the scatter plot
        thismodel = da1.sel(model=key).model.item()  # Label for legend
        thisplot = ax.scatter(da1.sel(model=key), da2.sel(model=key), marker=mymarkers[i], color=mycolors[i % 4],
                                  s=150, linewidths=1, label=thismodel)
    ax.scatter(da1_foc,da2_foc,marker='*',color='yellow',edgecolor='k',s=250,label=da1_foc.model.item())

    # Plot the line of best fit
    fitcorr = stats.pearsonr(da1.values, da2.values)
    yfit = np.poly1d(np.polyfit(da1.values, da2.values, 1))
    ax.plot(da1.values, yfit(da1.values), color='k', label='r={:.2f} (p={:.1e})'.format(fitcorr[0].item(),fitcorr[1].item()))
    handles, labels = ax.get_legend_handles_labels()
    
    #Labels
    ax.set_title('Climatological Bias compared to OMIP{}'.format(da1.OMIP.values), fontsize=font_title)
    ax.set_xlabel('{} ({}), [{}]'.format(da1.name,da1.units,da1.region), fontsize=font_label)
    ax.set_ylabel('{} ({}), [{}]'.format(da2.name,da2.units,da2.region), fontsize=font_label)
    ax.tick_params(axis='both', labelsize=font_tick)
    ax.grid(True)
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

    return handles,labels

def region_str(region):
    if (region[0]>180):
        x0 = str(360-region[0])+r'$^{\circ}$W'
    else:
        x0 = str(region[0])+r'$^{\circ}$E'
    if (region[1]>180):
        x1 = str(360-region[1])+r'$^{\circ}$W'
    else:
        x1 = str(region[1])+r'$^{\circ}$E'
    if (region[2]>0):
        y0 = str(region[2])+r'$^{\circ}$N'
    else:
        y0 = str(abs(region[2]))+r'$^{\circ}$S'
    if (region[3]>0):
        y1 = str(region[3])+r'$^{\circ}$N'
    else:
        y1 = str(abs(region[3]))+r'$^{\circ}$S'
    regstr = '{}-{}, {}-{}'.format(x0,x1,y0,y1)
    return regstr

def error_stats(da,region,fregion):
    varname = da.name
    ds_out = da.rename(varname+'_bias').to_dataset()
    regweights = np.cos(np.deg2rad(da.lat)) #.drop_attrs()
    da_freg = da.sel(lat=slice(fregion[2], fregion[3]), lon=slice(fregion[0], fregion[1]))
    fregweights = np.cos(np.deg2rad(da_freg.lat)) #.drop_attrs()
    ds_out[varname+'_bias'] = da_freg.weighted(fregweights).mean(['lat', 'lon'], skipna=True, keep_attrs=True)
    ds_out[varname+'_bias'] = ds_out[varname+'_bias'].assign_attrs({'description':'Focus region mean bias',
                                                                   'region':region_str(fregion)})
    ds_out[varname+'_rmse'] = ((da**2).weighted(regweights).mean(['lat', 'lon'], skipna=True))**0.5
    ds_out[varname+'_rmse'] = ds_out[varname+'_rmse'].assign_attrs({'description':'Plot region RMSE','units':da.units,'time_avg':da.time_avg,
                                                                   'region':region_str(region)})
    return ds_out

def plot_preproc(da,region,month):
    monstr = {1:'JAN',2:'FEB',3:'MAR',4:'APR',5:'MAY',6:'JUN',7:'JUL',8:'AUG',9:'SEP',10:'OCT',11:'NOV',12:'DEC'}
    da = da.sel(lat=slice(region[2], region[3]), lon=slice(region[0], region[1]))
    
    # Add 'month' coordinate from 'time'
    if 'time' in da.coords:
        da = da.assign_coords(month=da['time'].dt.month)
        if (month==13):
            # annual mean
            da = da.groupby('month').mean('time').mean('month').assign_attrs({'time_avg': 'Annual'})
            # da = da.mean('month').assign_attrs({'time_avg':'Annual'})
        else:
            # monthly mean
            da = da.groupby('month').mean('time').sel(month=month).assign_attrs({'time_avg': monstr[month]})
            # da = da.sel(month=month).assign_attrs({'time_avg':monstr[month]})
    elif 'month' in da.dims:
        # Already has monthly mean as a dimension
        if month == 13:
            da = da.mean('month').assign_attrs({'time_avg': 'Annual'})
        else:
            da = da.sel(month=month).assign_attrs({'time_avg': monstr[month]})

    else:
        raise ValueError("DataArray must have either 'time' coordinate or 'month' dimension for time averaging")

    return da

def get_units(var):
    unit_dict = {'thetao_zavg':r'$^{\circ}$C','thetao':r'$^{\circ}$C','so_zavg':'psu','so':'psu','sigma0':r'kg m$^{-3}$','sigma0_zavg':r'kg m$^{-3}$','mld':'m'}
    return unit_dict[var]

def SpatialPlot_climo_bias(ds_target, ds_model, ds_obs, var, region=[360-90, 360-0, 20, 80], focus_region=[360-48, 360-30, 38, 53], month=13, save=False, savedir='./'):
    '''
    SpatialPlot_climo_bias generates spatial plots of climatological bias and bias rank relative to a CMIP6 OMIP2 simulation library.

    Parameters:
        ds_target (xarray.Dataset) : Dataset containing the target model data to be analyzed. Must have dimensions: 'model', 'lat', 'lon', 'month'
        ds_model (xarray.Dataset) : Dataset containing the OMIP2 model data library. Must have dimensions: 'model', 'lat', 'lon', 'month'
        ds_obs (xarray.Dataset) : Dataset containing the benchmark observational data. Must have dimensions: 'lat', 'lon', 'month'
        var (str) : The variable to be analyzed. (Select from: ['thetao_zavg', 'so_zavg', 'sigma0_zavg', 'mld'])
        region (list) : The spatial domain for the plot in format [x0,x1,y0,y1] where x0,x1 span [0..360] and y0,y1 span [-90..90]
        focus_region (list) : The spatial subdomain computing regional bias, in format [x0,x1,y0,y1] where x0,x1 span [0..360] and y0,y1 span [-90..90]
        month (int) : The month of climatology to analyze. Set to 13 for annual-average.
        save (bool) : Save figures if True.

    Returns:
        ds_stats (xarray.Dataset) : Dataset containing area-weighted RMSE (over full plot domain) and regional bias (over focus region domain).
    '''
    da_target = plot_preproc(ds_target[var],region,month)
    da_model = plot_preproc(ds_model[var],region,month)
    da_obs = plot_preproc(ds_obs[var],region,month)
    
    # Add target to models
    focus_model = ds_target.model.item()
    da_model = xr.concat([da_model, da_target], dim='model')

    # Compute Bias
    da_model = (da_model - da_obs).assign_attrs({'units':get_units(var),'time_avg':da_obs.time_avg})
    ds_stats = error_stats(da_model,region,focus_region)

    # Plot
    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(20, 10), layout='constrained')
    gs = GridSpec(1, 2, figure=fig)

    axs11 = fig.add_subplot(gs[0, 0], projection=proj)
    SpatialBias_panel(da_model, ds_stats, focus_region=focus_region, focus_model=focus_model, ax=axs11)
        
    axs12 = fig.add_subplot(gs[0, 1], projection=proj)
    SpatialRank_panel(da_model, focus_region=focus_region, focus_model=focus_model, ax=axs12)

    # Save Plots
    if save:
        plotname = savedir+'/SpatialPlot_climo_bias.{}.{}.{}.png'.format(focus_model.replace(" ", "_"),var,da_model.time_avg)
        fig.savefig(plotname)

    return ds_stats

def ScatterPlot_Error(ds_x, var_x, ds_y, var_y, focus_model, save=False, savedir='./'):
    '''
    ScatterPlot_Error generates a scatter plot relating error statistics generated through calls to SpatialPlot_climo_bias().

    Parameters:
        ds_x (xarray.Dataset) : Dataset containing the error statistics to be analyzed. Must have dimensions: 'model'
        ds_y (xarray.Dataset) : Dataset containing the error statistics to be analyzed. Must have dimensions: 'model'
        focus_model (str) : Name of focus model 
    '''

    fig = plt.figure(figsize=(10, 5), layout='constrained')
    gs = GridSpec(1, 2, figure=fig)
    
    # Add Spatial Bias Plots
    axs11 = fig.add_subplot(gs[0, 0], label='var1 bias')
    han,lab = Scatter_panel(ds_x[var_x], ds_y[var_y], focus_model, axs11)
    axs12 = fig.add_subplot(gs[0, 1], label='legend')
    axs12.legend(han, lab, ncols=1, loc='center', fontsize=12, markerscale=0.75)  # Add legend to the dedicated axes
    axs12.axis('off')

    # Save Plots
    if save:
        plotname = savedir+'/ScatterPlot_Error.{}.{}_{}.{}_{}.png'.format(focus_model.replace(" ", "_"),var_x,ds_x[var_x].time_avg,var_y,ds_y[var_y].time_avg)
        fig.savefig(plotname)

    return 



#### PART 3 CALCULATION


def compute_wmt(ds, calc_type, density, dsigma, dclasses=np.empty(0), verbose=True):
    """
    Computes water mass transformation in the subpolar North Atlantic.

    Parameters:
    ds: xarray.Dataset
         Dataset with fields necessary for WMT calculation (TEMP, SALT, WFO/VSF, HFDS)
    calc_type: string
         Keyword with the type of calculation to perform
         - WMT: line plots of WMT vs sigma
         - MAPS: transformation in sigma decomposed also by location
         - DENS: ??  
         - DS:  ??
    density: string
         Keyword with sigma variable to use (i.e., sigma2 or sigma0)
         sigma2 implemented in the POD for comparison to AMOC(sigma2)
    dsigma: float
         width of sigma bins
    dclasses: numpy array, optional
         array with water mass values to use for WMT calculations
    verbose: boolean, optional

    Returns:
    xarray.Dataset
        Dataset with water mass transformation in SPNA regions. 
    """
    if calc_type == 'WMT':
        ds_wmt = calc_wmt(ds, density, dclasses, dsigma)
    if calc_type == 'MAPS':
        ds_wmt = calc_maps(ds, density, dclasses, dsigma)
    if calc_type == 'DENS':
        ds_wmt = calc_dens(ds, density, dclasses)
    if calc_type == 'DS':
        ds_wmt = make_ds(ds)

    return ds_wmt


def make_spna_masks(ds):
    """
    Makes masks of subpolar North Atlantic subregions 

    Parameters:
    ds: xarray.Dataset
        Dataset that includes lon and lat needed to determine region bounds

    Returns:
    xarray.Dataset
        Dataset with region masks applied

    """
    dbasins = spna_masks.spna_masks()
    lon_raw = ds['lon']  # TO-DO may need un-hardcoding
    lat = ds['lat']  #TO-DO may need un-hardcoding
    lon = xr.where(lon_raw < -180, lon_raw + 360.0, lon_raw)
    #having the wraparound point in atlantic breaks masks for some reason
    lon = xr.where(lon > 180, lon-360.0, lon)
    #fix for 1-d lat/lon
    if len(lon.dims)<2:
        lonlon, latlat =np.meshgrid(lon,lat)
        reg_mask = dbasins.mask(lonlon, latlat)
    else:
        reg_mask = dbasins.mask(lon, lat)
    lab = ds.where(reg_mask == 1)
    spg_sw = ds.where(reg_mask == 2)
    spg_se = ds.where(reg_mask == 3)
    irm = ds.where(reg_mask == 4)
    nor = ds.where(reg_mask == 5)
    arc = ds.where(reg_mask == 6)
    spna = ds.where(reg_mask < 6)
    ds_mask = xr.concat([lab,spg_sw,spg_se,irm,nor,arc,spna], 'region')
    regions = ['Labrador Sea','Subpolar Gyre SW', 'Subpolar Gyre SE', 'Irminger-Icelandic Sea',
                   'Nordic Seas', 'Arctic Sea', 'Subpolar North Atlantic']
    ds_mask = ds_mask.assign_coords({'region':regions})

    return ds_mask

def wmt_preproc(ds):
    """
    Preprocessing of Dataset for xwmt functions

    Parameters: 
    ds: xarray.Dataset
        Dataset that includes all fields needed for WMT calculation

    Returns:
    xarray.Dataset   
        Dataset with fields ready for xwmt processing

    """
    # Surface freshwater (mass) flux for xwmt. CESM provides a virtual salt flux
    # (vsf); most other models provide wfo directly. Prefer wfo when present;
    # otherwise convert vsf. (When both exist, the upstream loader passes only
    # wfo, but the same preference is enforced here too.)
    if 'wfo' in ds:
        # CMIP `wfo` (water_flux_into_sea_water) is positive INTO the ocean,
        # which is exactly the convention xwmt expects (it forms the salt forcing
        # as -wfo*sos), so it is used as-is. This matches the sign of the vsf->wfo
        # conversion below (vsf>0 = salt in -> wfo<0 = freshwater out). No -1
        # flip: that would make the native-wfo path inconsistent with the vsf path.
        pass
    elif 'vsf' in ds:
        # Virtual salt flux -> freshwater (mass) flux. xwmt expects `wfo` to be a
        # freshwater mass flux (kg m-2 s-1), NOT a salt flux, so vsf cannot be
        # assigned to wfo directly. FSU-style conversion from HR-LR_OSNAP_wmt
        # (sub2sub_HR/wmt.py::_op_sf_to_fwf_from_sos):
        #     wfo = vsf * (-1 / sos)
        # No 1000 factor: xwmt pairs wfo with sos as a -wfo*sos salt forcing, so
        # vsf [kg m-2 s-1] divided by sos [psu] lands in the same convention.
        if 'sos' not in ds:
            raise KeyError("wmt_preproc: 'sos' is required to convert 'vsf' -> 'wfo'")
        ds['wfo'] = ds['vsf'] * (-1.0 / ds['sos'])
    else:
        raise KeyError("wmt_preproc: need either 'wfo' or 'vsf' to build the freshwater flux for WMT")
    ds['wet'] = xr.where(~np.isnan(ds.tos), 1, 0)
    ds['sfdsi'] = xr.zeros_like(ds['hfds']).rename('sfdsi')
    return ds


def make_ds(ds):
    """
    Applies all preprocessing needed for regions and xwmt conventions

    Parameters:
    ds: xarray.Dataset
        Dataset that includes fields needed for calculation

    Returns:
    xarray.Dataset

    """
    ds_preproc = wmt_preproc(ds)
    ds_xwmt = make_spna_masks(ds_preproc)
    return ds_xwmt

def calc_wmt(ds, density, dclasses, dsigma):
    """
    calculation of WMT for entire regions using xwmt 

    Parameters:
    ds: xarray.Dataset
        Preprocessed dataset with all fields needed for WMT
    density: string
        Keyword that describes which xwmt function will be called
    dclasses: numpy array
        User-defined array of water mass classes        
    dsigma: float
        Width of sigma bins

    Returns: 
    xarray.Dataset
        Dataset with WMT as a function of sigma only
    """

    #xwmt_ds = make_ds(ds, region=region)
    ds_preproc = wmt_preproc(ds)
    ds_mask = make_spna_masks(ds_preproc)
    xwmt_init = xwmt.swmt(ds_mask)
    if (density == 'sigma0') and (len(dclasses)==0):
        bins = np.arange(24.5, 29.6, dsigma)
    elif (density == 'sigma2') and (len(dclasses)==0):
        bins = np.arange(33.1, 38.6, dsigma)
    else: 
        bins = dclasses
    ds_wmt = xwmt_init.G(density, bins=bins)
    ds_wmt = ds_wmt.to_dataset(name='wmt')
    ds_wmt_decomp = xwmt_init.G(density, bins=bins, group_tend=False)
    print(ds_wmt_decomp)
    #ds_wmt['heat'] = ds_wmt_decomp['heat']
    #ds_wmt['freshwater'] = ds_wmt_decomp['freshwater']
    return ds_wmt

def calc_maps(ds, density, dclasses, dsigma):
    """
    calculation of WMT in each grid cell as a function of sigma using xwmt 

    Parameters:
    ds: xarray.Dataset
        Preprocessed dataset with all fields needed for WMT
    density: string
        Keyword that describes which xwmt function will be called
    dclasses: numpy array
        User-defined values of water masses for maps  
    dsigma: float
        Width of sigma bins  

    Returns: 
    xarray.Dataset
        Dataset with WMT as a function of sigma and location
    """

    print(dclasses, len(dclasses))
    ds_preproc = wmt_preproc(ds)
    ds_mask = make_spna_masks(ds_preproc)
    ds_mask.time.attrs['calendar_type'] = 'noleap'
    if (density == 'sigma0') and (len(dclasses)==0):
        vals = np.arange(26, 28.6, dsigma)
    elif (density == 'sigma2') and (len(dclasses)==0):
        vals = np.arange(35.5,38,dsigma)
    else: 
        vals = dclasses
    xwmt_init = xwmt.swmt(ds_mask.sel(region='Subpolar North Atlantic'))
   
    temp_array = []
    for dd in vals:
        print(dd) 
        trans1 =  xwmt_init.F(density, bins = np.array([dd-dsigma/2, dd+dsigma/2]), group_tend=True)
        temp_array.append(trans1) 
        
    wmt_maps_decomp = xr.concat(temp_array, dim='sigma2')
#    wmt_maps_decomp['total'] = wmt_maps
    return wmt_maps_decomp

def calc_dens(ds, density):
    """
    calculation of density using xwmt 

    Parameters:
    ds: xarray.Dataset
        Preprocessed dataset with all fields needed for WMT
    density: string
        Keyword that describes which xwmt function will be called
         

    Returns: 
    xarray.Dataset
        Dataset with density
    """

    ds_preproc = wmt_preproc(ds)
    ds_mask = make_spna_masks(ds_preproc)
    xwmt_init = xwmt.swmt(ds_mask)
    xwmt_den = xwmt_init.get_density(density)
    xwmt_den_ar = xwmt_den[2].isel(lev_outer=0)
    xwmt_den_ds = xwmt_den_ar.to_dataset(name=density)
    return xwmt_den_ds


#### PART 3 PLOTTING

def wmt_plot_byregion(ds_benchmark, ds_model, sigma_classes, save=False, savedir='./'):
    """
    POD plot of WMT lines by region in model versus observational benchmarks

    Parameters:
    ds_benchmark: xarray.Dataset
        dataset with observational benchmarks from Low et al. 
    ds_model: xarray.Dataset
        dataset with wmt in model simulation calculated using xwmt
    sigma_classes: list
        lower bound of sigma classes of interest for extra analysis 
    save: boolean, optional
        save figure or not for POD 
    savedir: string, optional
        location where figure is saved


    """
    #get time mean wmt from obs benchmarks
    wmt_benchmark =  ds_benchmark['wmt']

    #mean, min, max of WMT in each region needed for plot
    wmt_benchmark_mean = wmt_benchmark.mean('benchmark')
    wmt_benchmark_min = wmt_benchmark.min('benchmark')
    wmt_benchmark_max = wmt_benchmark.max('benchmark')

    #get time mean model  wmt
    wmt_model = ds_model['wmt']/1e6

    
    #add wmt plot
    fig = plt.figure(figsize=(6,10), layout='constrained')

    selected_regions = ['Subpolar North Atlantic', 'Nordic Sea','Labrador Sea', 'Subpolar Gyre SW', 'Irminger-Icelandic Sea', 'Subpolar Gyre SE']
    region_dim = wmt_benchmark.region
    
    numreg = len(selected_regions)
    gs = GridSpec(numreg, 1, figure=fig)

    #add WMT by region

    minwmt = min(wmt_benchmark_min.min(['region', 'sigma2']), wmt_model.min(['region', 'sigma2'])).values.item()
    maxwmt = max(wmt_benchmark_max.max(['region', 'sigma2']), wmt_model.max(['region', 'sigma2'])).values.item()
    buffer = (maxwmt-minwmt)/20
    minwmt = minwmt-buffer
    maxwmt = maxwmt+buffer

    minsig = min(wmt_benchmark.sigma2.min(), wmt_model.sigma2.min()).values.item()
    maxsig = max(wmt_benchmark.sigma2.max(), wmt_model.sigma2.max()).values.item()
    
    for rr in range(numreg):
        reg = selected_regions[rr]
        bool_reg = region_dim == reg
        
        ax = fig.add_subplot(gs[rr,0])
        ax.plot(wmt_benchmark_mean.sigma2, wmt_benchmark_mean.sel(region=bool_reg).squeeze(), color='black', label='benchmarks')
        ax.fill_between(wmt_benchmark_mean.sigma2, wmt_benchmark_min.sel(region=bool_reg).squeeze(), wmt_benchmark_max.sel(region=bool_reg).squeeze(),  alpha=0.5, color='gray')
        ax.plot(wmt_model.sigma2, wmt_model.sel(region=bool_reg).squeeze(), color='red', label='model')
        if rr==0: ax.legend(loc='upper right')
        ax.text(0.05, 0.8, reg, transform=ax.transAxes)
        
        ax.set_ylim(minwmt, maxwmt)
        ax.set_xlim(minsig, maxsig)
        
        ax.set_ylabel('WMT (Sv)')
        ax.axhline(0, color='gray', alpha=0.8, linewidth=1)
        [ax.axhline(ss, color='green', alpha=0.8, linewidth=2) for ss in sigma_classes]
        ax.grid(color='gray', linewidth=1, linestyle='dashed', alpha=0.5)
        #[ax.axvline(vv, color='gray', alpha=0.8, linewidth=1) for vv in np.arange(int(np.floor(minsig))+1,int(np.ceil(maxsig)),1)]
        if rr == numreg-1: 
            ax.set_xlabel('$\\sigma_{2}$ (kg m$^{-3}$)')
   
    #Save Plots
    if save:
        plotname = savedir+'/wmt_lineplot_byregion.png'
        plt.savefig(plotname)

    return



def wmt_amoc_plot(ds_wmt_benchmarks, ds_wmt_model, ds_moc, lat_target=45,
                  region_name='Subpolar North Atlantic', save=False, savedir='./'):
    """
    Three-curve comparison vs sigma2 in the SPNA: obs WMT (mean + spread),
    model WMT, and model AMOC at a target latitude. Mass conservation
    implies AMOC(lat_target, sigma2) ~ -WMT integrated over SPNA(lat>lat_target)
    in steady state, so these curves should approximately compensate.

    Parameters
    ----------
    ds_wmt_benchmarks : xarray.Dataset
        Observation-based WMT benchmarks (dims: benchmark, region, sigma2;
        var: wmt). Loaded from obs_wmt_sigma2_*.nc.
    ds_wmt_model : xarray.Dataset
        Model WMT from POD_utils.compute_wmt (dims: region, sigma2; var: wmt,
        units m^3/s -- divided by 1e6 here to get Sv).
    ds_moc : xarray.Dataset
        Model AMOC from POD_utils.calculate_moc (dims: region, sigma, lat,
        time; var: MOC in Sv). region=1 is Atlantic+Arctic.
    lat_target : float
        Latitude (deg N) at which to extract the AMOC curve. Default 45.
    region_name : str
        WMT region label to select for both obs and model. Default
        'Subpolar North Atlantic' (the combined SPNA total).
    save, savedir : bool, str
        If save=True, write the figure to savedir/wmt_amoc{lat_target}.png.
    """
    # Obs WMT for the selected region: mean line + min/max benchmark spread
    region_dim = ds_wmt_benchmarks.region
    bool_reg = region_dim == region_name
    wmt_obs_mean = ds_wmt_benchmarks['wmt'].mean('benchmark').sel(region=bool_reg).squeeze()
    wmt_obs_min = ds_wmt_benchmarks['wmt'].min('benchmark').sel(region=bool_reg).squeeze()
    wmt_obs_max = ds_wmt_benchmarks['wmt'].max('benchmark').sel(region=bool_reg).squeeze()

    # Model WMT for the same region; convert m^3/s -> Sv
    wmt_model = ds_wmt_model['wmt'].sel(region=bool_reg).squeeze() / 1e6

    # Model AMOC at lat_target, Atlantic+Arctic basin, time mean.
    # .sel(method='nearest') may snap to 44.5 or 45.5 on the 1-deg grid.
    moc_at_lat = ds_moc['MOC'].isel(region=1).sel(lat=lat_target, method='nearest').mean('time')

    # Streamfunction type set by calculate_moc: 'Eulerian' (uo/vo) or 'residual'
    # (umo/vmo). Fall back to a generic label if the attr is missing.
    sf_type = ds_moc.attrs.get('streamfunction_type',
                               ds_moc['MOC'].attrs.get('streamfunction_type'))
    sf_label = f'{sf_type} streamfunction' if sf_type else 'AMOC'

    fig, ax = plt.subplots(figsize=(8, 5))

    # Obs benchmark spread (gray fill) + mean (black line) -- matches wmt_plot_byregion style
    ax.fill_between(wmt_obs_mean.sigma2, wmt_obs_min, wmt_obs_max,
                    alpha=0.4, color='gray', label='Obs WMT benchmarks range')
    ax.plot(wmt_obs_mean.sigma2, wmt_obs_mean, color='black', label='Obs WMT (mean)')

    # Model WMT (red, same color as wmt_plot_byregion)
    ax.plot(wmt_model.sigma2, wmt_model, color='red', label='Model WMT')

    # Model AMOC at target latitude (blue). ds_moc's coord is named 'sigma' (not
    # 'sigma2'); the values are sigma2 by construction in calculate_moc.
    ax.plot(moc_at_lat.sigma, moc_at_lat, color='blue',
            label=f'Model {sf_label} at {lat_target}$^\\circ$N')

    # Limit x-axis to physical sigma2 range (the ds_moc 'sigma' coord includes
    # a 0 anchor from sigma2_grid_96L that would otherwise compress the view).
    sig_lo = float(wmt_obs_mean.sigma2.min())
    sig_hi = float(wmt_obs_mean.sigma2.max())
    ax.set_xlim(sig_lo, sig_hi)

    ax.axhline(0, color='gray', alpha=0.6, linewidth=1)
    ax.grid(color='gray', linewidth=1, linestyle='dashed', alpha=0.5)
    ax.set_xlabel(r'$\sigma_2$ (kg/m$^3$)')
    ax.set_ylabel('Volume flux (Sv)')
    ax.set_title(f'{region_name}: WMT vs {sf_label} at {lat_target}$^\\circ$N')
    ax.legend(loc='best')
    plt.tight_layout()

    if save:
        plotname = f'{savedir}/wmt_amoc{lat_target}.png'
        plt.savefig(plotname)

    return



def wmt_plot_maps(ds_benchmark, ds_model, dimnames, sigma_classes, save=False, savedir='./'):
    """
    POD plot of WMT lines by region in model versus observational benchmarks

    Parameters:
    ds_benchmark: xarray.Dataset
        dataset with observation-based benchmarks from Low et al. 
    ds_model: xarray.Dataset
        dataset with sigma in model simulation calculated using xwmt
    dimnames: list
        list that includes dimension names for model output
    sigma_classes: list
        lower bound of sigma classes of interest for maps 
    save: boolean, optional
        save figure or not for POD 
    savedir: string, optional
        location where figure is saved


    """

    nsigma = len(sigma_classes)
    gs=GridSpec(nsigma,2)

    time_coord = dimnames[0]
    lon_coord = dimnames[1]
    lat_coord = dimnames[2]

    #correct lon for North Atlantic
    model_lon = ds_model['lon'] #ok to hard-code 'lon' and 'lat' here b/c this is output from xwmt
    if model_lon.max()>345:
        print('correcting lon for dateline')
        # Build a wrapped copy with xr.where instead of mutating model_lon.values
        # in place -- the regridded lon array is read-only, so the in-place write
        # raised "assignment destination is read-only".
        model_lon = xr.where(model_lon>180, model_lon-360.0, model_lon)
        ds_model = ds_model.assign_coords({'lon':model_lon})
    model_lat = ds_model['lat']


    #making mpl colormap from ncl colortable
    cmap_name='nrl_sirkes'
    colortab=pd.read_csv('https://www.ncl.ucar.edu/Document/Graphics/ColorTables/Files/'+cmap_name+'.rgb',sep='\\s+')
    cmap = LinearSegmentedColormap.from_list(cmap_name, colortab.values/255, N=101)
    cmap.set_bad(color='white')

    # Shared, symmetric color limits derived from BOTH fields (99th pct of
    # |trans.|), so the model and benchmark panels sit on the same data-driven
    # scale instead of a range hard-tuned to the obs benchmark.
    _mod_plot = ds_model.sel(sigma2=sigma_classes, method='nearest').mean(time_coord) / 1e6
    _ben_plot = ds_benchmark['wmt'].sel(sigma2=sigma_classes, method='nearest').mean('benchmark') / 1e6
    _mod_abs = np.abs(_mod_plot.values)
    vlim = float(np.nanmax([
        np.nanpercentile(_mod_abs, 99),
        np.nanpercentile(np.abs(_ben_plot.values), 99),
    ]))
    if not np.isfinite(vlim) or vlim == 0:
        vlim = 2e-11

    sigma2_mask_ben = ds_benchmark['wmt_freq'].mean('benchmark')
    # Model 'outcrop frequency': fraction of timesteps a cell actively transforms
    # this class. Threshold is relative to the field's own scale -- the old 1e-11
    # was ~1e6x too small for these ~1e-5 values, so it flagged every cell (incl.
    # NaNs outside the SPNA) as outcropping.
    _raw_scale = float(np.nanmax(_mod_abs)) * 1e6
    _mod_thresh = 1e-3 * _raw_scale if (np.isfinite(_raw_scale) and _raw_scale > 0) else 0.0
    sigma2_mask_mod = xr.where(np.abs(ds_model) > _mod_thresh, 1, 0)

    f=plt.figure(figsize=(16,2.25*nsigma))
    #loop through density classes of interest
    for ii,ss in enumerate(sigma_classes):
        #obs outcrop frequency
        obs_outcrop = sigma2_mask_ben.sel(sigma2=ss, method='nearest')
        #model outcrop frequency
        mod_outcrop = sigma2_mask_mod.sel(sigma2=ss, method='nearest').mean(time_coord) 

        #model plots
        ax=plt.subplot(gs[ii,0],projection = ccrs.PlateCarree())    
        cs2=plt.pcolormesh(model_lon, model_lat, ds_model.sel(sigma2=ss,method='nearest').mean(time_coord)/1e6, cmap=cmap, vmin=-vlim, vmax=vlim, transform=ccrs.PlateCarree())
    #cs=plt.contour(model_lon, model_lat, mod_outcrop, np.arange(0.1,0.31, 0.1), cmap=plt.cm.viridis, transform=ccrs.PlateCarree())
        plt.title('Model: '+'$\\sigma_{2}$='+str(ss)[0:4]+'-'+str(ss+0.1)[0:4]+' kg/m$^3$')
        plt.colorbar(cs2, label='trans. (Sv/m$^2$)')
        ax.coastlines()
        ax.set_extent([-80,30,45,80], crs=ccrs.PlateCarree())
    
        #benchmark plots
        ax=plt.subplot(gs[ii,1],projection = ccrs.PlateCarree())
        cs2=plt.pcolormesh(ds_benchmark.lon, ds_benchmark.lat, ds_benchmark['wmt'].sel(sigma2=ss,method='nearest').mean('benchmark')/1e6, cmap=cmap, vmin=-vlim, vmax=vlim, transform=ccrs.PlateCarree())
        cs=plt.contour(obs_outcrop.lon, obs_outcrop.lat, obs_outcrop, np.arange(0.1,0.31, 0.1),  cmap=plt.cm.viridis, transform=ccrs.PlateCarree())
        ax.set_extent([-80,30,45,80], crs=ccrs.PlateCarree()) 
        ax.coastlines()
        plt.title('Obs: '+'$\\sigma_{2}$='+str(ss)[0:4]+'-'+str(ss+0.1)[0:4]+' kg/m$^3$')
        plt.colorbar(cs2, label='trans. (Sv/m$^2$)')
    
    plt.tight_layout()

    #Save Plots
    if save:
        plotname = savedir+'/wmt_maps_selectclasses.png'
        plt.savefig(plotname)

    return

##### PART 4: AMOC #####
# Adapted from Sub2Sub's moc_funcs.py (Yeager, Maroon), at
# /glade/work/emaroon/Sub2Sub/sub2sub/moc_funcs.py. Computes the meridional
# overturning circulation streamfunction in sigma2 (potential density referenced
# to 2000 dbar) coordinates. Output dataset has MOC(region, sigma, lat, time)
# in Sverdrups; region indices are 0=Global, 1=Atlantic+Arctic, 2=IndoPac+SO.
#
# Public entry point: calculate_moc(ds_t, ds_u, ds_v, use_currents=False)
# - ds_t: tracer dataset, must contain thetao, so, lev, lev_bnds, time, time_bnds, lon, lat
# - ds_u: u-component dataset, must contain either uo (use_currents=True) or umo (False)
# - ds_v: v-component dataset, must contain either vo or vmo
# - use_currents: True uses uo/vo (velocity, converted to mass flux via geometry);
#                 False uses umo/vmo (mass transport, divided by reference density)

# Sigma-coord name (used internally for xhistogram bin naming) and Earth radius
# in meters (for great-circle distances).
_AMOC_SIGNAME = 'sigma'
_AMOC_EARTH_R = 6371e3


def which_grid(tlon, ulon, vlon):
    """
    Identify Arakawa grid type from U, V, and T-point longitudes.
    Returns a single letter: 'a' (collocated), 'b' (U==V, offset from T),
    'c' (V==T, U offset), 'd' (U==T, V offset), or 'other'.
    """
    if (ulon == vlon).all() & (ulon != tlon).any(): grid = 'b'
    elif (ulon != vlon).any() & (vlon == tlon).any(): grid = 'c'
    elif (ulon != vlon).any() & (ulon == tlon).all(): grid = 'd'
    elif (ulon == vlon).all() & (ulon == tlon).all(): grid = 'a'
    else: grid = 'other'
    return grid


def sigma2_grid_96L():
    """
    Define a 156-point sigma2 grid for MOC(sigma2) binning. Returns midpoints
    and edge points as DataArrays suitable for xhistogram.
    Despite the function name (kept from Sub2Sub), the grid has 156 layers
    spanning sigma2 = 26 to 38.
    """
    tmp1 = np.arange(26, 35, 0.2)
    tmp2 = np.arange(35, 36, 0.1)
    tmp3 = np.arange(36, 38.05, 0.05)
    sig2 = np.concatenate((tmp1, tmp2, tmp3))
    sigma_mid = xr.DataArray(
        sig2, coords={_AMOC_SIGNAME: sig2},
        attrs={'long_name': 'Sigma2 at middle of layer', 'units': 'kg/m^3'},
    )
    sigma_edge = (sigma_mid + sigma_mid.shift(sigma=1)) / 2.0
    sigma_edge[0] = 0.0
    sigma_edge = np.append(sigma_edge.values, [50.0])
    sigma_edge = xr.DataArray(
        sigma_edge, coords={_AMOC_SIGNAME: sigma_edge},
        attrs={'long_name': 'Sigma2 at edges of layer', 'units': 'kg/m^3'},
    )
    return sigma_mid, sigma_edge


def _amoc_fluxdiv_B(uflux, vflux):
    """B-grid horizontal flux divergence. Assumes uflux=U*DY*DZ, vflux=V*DX*DZ."""
    UTE = 0.5 * (uflux + uflux.shift(y=1))
    UTW = UTE.roll(x=1, roll_coords=False)
    VTN = 0.5 * (vflux + vflux.roll(x=1, roll_coords=False))
    VTS = VTN.shift(y=1)
    return UTE - UTW + VTN - VTS


def _amoc_fluxdiv_C(uflux, vflux):
    """C-grid horizontal flux divergence."""
    UTE = uflux
    UTW = UTE.roll(x=1, roll_coords=False)
    VTN = vflux
    VTS = VTN.shift(y=1)
    return UTE - UTW + VTN - VTS


def wflux_div(grid, uflux, vflux, densdim, densedges):
    """
    Vertical volume flux in density-space at T-point, derived from horizontal
    u/v fluxes by computing convergence (= -divergence) and integrating from
    the ocean bottom upward.

    Parameters
    ----------
    grid : str
        Arakawa grid type ('a', 'b', or 'c').
    uflux, vflux : xarray.DataArray
        Horizontal volume fluxes in density coordinates (m^3/s).
    densdim : str
        Name of the density dimension to integrate over.
    densedges : xarray.DataArray
        Density-layer edge values.

    Returns
    -------
    xarray.DataArray
        Vertical volume flux (m^3/s) at T-point, on the same density edges.
    """
    # Convergence on the chosen grid type
    if grid == 'b':
        dwflux = -_amoc_fluxdiv_B(uflux, vflux)
    else:
        dwflux = -_amoc_fluxdiv_C(uflux, vflux)

    # Bottom-up vertical (density) integral to recover W
    kwargs = {densdim: slice(None, None, -1)}
    wflux = dwflux.sel(kwargs).cumsum(densdim).sel(kwargs)

    kwargs = {densdim: slice(0, -1)}
    wflux['sigma'] = densedges.isel(kwargs)
    return wflux


def latitude_grid_1deg():
    """1° latitude grid for MOC output. Returns (mid, edge) DataArrays."""
    midvals = np.arange(-89.5, 90.5, 1)
    edgevals = np.arange(-90, 91, 1)
    lat_mid = xr.DataArray(
        midvals, coords={'lat': midvals},
        attrs={'long_name': 'latitude', 'units': 'degrees_north'}, name='lat',
    )
    lat_edge = xr.DataArray(
        edgevals, coords={'lat': edgevals},
        attrs={'long_name': 'latitude', 'units': 'degrees_north'}, name='lat_edge',
    )
    return lat_mid, lat_edge


# Cached at import time — used as the default output grid for calculate_moc
_AMOC_LAT_MID, _AMOC_LAT_EDGE = latitude_grid_1deg()


def basinmask(da, lon_name='lon', lat_name='lat'):
    """
    Build ocean-basin masks (Global, Atlantic+Arctic, Indo-Pacific+SO) from
    cmip_basins.basins. Returns a region-dimensioned DataArray of 0/1 masks.
    """
    grid = xr.Dataset()
    if len(da[lon_name].shape) == 1:
        xx, yy = np.meshgrid(da[lon_name], da[lat_name])
        grid.coords['lon'] = xr.DataArray(xx, dims=['y', 'x'], name='lon')
        grid.coords['lat'] = xr.DataArray(yy, dims=['y', 'x'], name='lat')
    else:
        grid.coords['lon'] = da[lon_name]
        grid.coords['lat'] = da[lat_name]

    codes = cmip_basins.basins.generate_basin_codes(
        grid, lon='lon', lat='lat', persian=False, style='cmip6',
    )

    global_codes = np.arange(11).tolist()
    atl_codes = [2, 4, 6, 7, 8, 9]
    indopac_codes = [1, 3, 5, 10]
    atl_mask = xr.concat([xr.where(codes == aa, 1, 0) for aa in atl_codes], dim='dummy').sum('dummy')
    indopac_mask = xr.concat([xr.where(codes == ii, 1, 0) for ii in indopac_codes], dim='dummy').sum('dummy')
    global_mask = xr.concat([xr.where(codes == ii, 1, 0) for ii in global_codes], dim='dummy').sum('dummy')

    mask = xr.concat([global_mask, atl_mask, indopac_mask], dim='region')
    mask.attrs['legend'] = {0: 'Global', 1: 'Atlantic+Arctic', 2: 'IndoPac+SO'}
    return mask


def wflux_zonal_sum(wflux, tlat, regionmask, lat, lat_bin_name='lat_t_bin'):
    """
    Zonally integrate w-flux per basin using xhistogram binning by latitude.

    Parameters
    ----------
    wflux : xarray.DataArray
        Vertical volume flux (m^3/s).
    tlat : xarray.DataArray
        T-grid latitude (2D).
    regionmask : xarray.DataArray
        Per-region 0/1 masks (region dim).
    lat : xarray.DataArray
        Target latitude bins (mid-points).
    lat_bin_name : str
        Name xhistogram assigns to the binned latitude axis.

    Returns
    -------
    xarray.DataArray
        Zonally-integrated wflux per region on the target lat grid.
    """
    wgts = (wflux * regionmask).astype('float32')
    xr_out = histogram(tlat, bins=[lat.data], weights=wgts, dim=['y', 'x'], density=False)

    # Zero at southern edge to prepare for meridional integral
    xr_out[{lat_bin_name: 0}] = 0
    xr_out = xr_out.rename({lat_bin_name: lat.name})
    xr_out = xr_out.assign_coords({_AMOC_SIGNAME: wflux[_AMOC_SIGNAME]})
    xr_out[lat.name] = lat[1:]
    return xr_out


def compute_MOC(wflux, tlat, regionmask, lat, lat_bin_name='lat_bin'):
    """
    W-method MOC: zonally sum wflux per basin (xhistogram by latitude), then
    cumulatively sum meridionally from south to north. Result in Sverdrups.
    """
    zonsum = wflux_zonal_sum(wflux, tlat, regionmask, lat, lat_bin_name=lat_bin_name)
    moc = zonsum.cumsum(dim=lat.name) / 1.0e6
    moc = moc.assign_attrs({'long_name': 'Meridional Overturning Circulation', 'units': 'Sv'})
    moc.name = 'MOC'
    return moc


def _amoc_hav(x):
    """Haversine helper: sin(x/2)^2."""
    return np.sin(x / 2) ** 2


def _amoc_archav(x):
    """Inverse-haversine helper: 2*arcsin(sqrt(x))."""
    return 2 * np.arcsin(np.sqrt(x))


def great_circ_dist2(lat, lon, dimname):
    """
    Haversine great-circle distance between adjacent grid points in dimname
    ('x' uses roll, 'y' uses shift since the y boundaries don't wrap).
    Returns distance in meters; NaNs out where lat is unphysical (>1e30).
    """
    if dimname == 'x':
        dlat = (lat - lat.roll(x=1)) * np.pi / 180
        dlon = (lon - lon.roll({dimname: 1})) * np.pi / 180
        inside = _amoc_archav(
            _amoc_hav(dlat)
            + (1 - _amoc_hav(dlat) - _amoc_hav((lat.roll(x=1) + lat) * np.pi / 180))
            * _amoc_hav(dlon)
        )
    elif dimname == 'y':
        dlat = (lat - lat.shift({dimname: 1})) * np.pi / 180
        dlon = (lon - lon.shift({dimname: 1})) * np.pi / 180
        inside = _amoc_archav(
            _amoc_hav(dlat)
            + (1 - _amoc_hav(dlat) - _amoc_hav((lat.shift({dimname: 1}) + lat) * np.pi / 180))
            * _amoc_hav(dlon)
        )
        # First column has no left neighbor — copy from second column to avoid NaN
        inside.values[:, 0] = inside.values[:, 1]
    else:
        raise ValueError(f"great_circ_dist2: dimname must be 'x' or 'y', got {dimname!r}")

    dsig = inside.where(lat < 1e30)
    return _AMOC_EARTH_R * dsig


def calculate_moc(ds_t, ds_u, ds_v, use_currents=False):
    """
    Compute the meridional overturning circulation streamfunction in sigma2
    coordinates from three input datasets.

    Parameters
    ----------
    ds_t : xarray.Dataset
        Tracer dataset containing thetao, so, lev, lev_bnds, time, time_bnds,
        lon, lat. Must already be on y/x dim names with lev in meters
        (run POD_utils.preprocess_coords first).
    ds_u, ds_v : xarray.Dataset
        Velocity datasets. If use_currents=True, each must contain uo / vo
        (m/s). If use_currents=False, must contain umo / vmo (kg/s).
    use_currents : bool
        True  -> velocity path (uo, vo + grid geometry -> volume flux)
        False -> mass-transport path (umo, vmo / 1028 kg/m^3 -> volume flux)

    Returns
    -------
    xarray.Dataset
        Contains 'MOC' (region, sigma, lat, time) in Sverdrups, plus
        'time_bnds'. region=0 Global, 1=Atlantic+Arctic, 2=IndoPac+SO.
    """
    # Drop conflicting x/y scalar coords if present (preprocess_coords leftovers)
    for coord in ('x', 'y'):
        if coord in ds_t.coords: ds_t = ds_t.drop_vars(coord)
        if coord in ds_u.coords: ds_u = ds_u.drop_vars(coord)
        if coord in ds_v.coords: ds_v = ds_v.drop_vars(coord)

    sigma_mid, sigma_edge = sigma2_grid_96L()

    tlon, tlat = ds_t['lon'], ds_t['lat']
    ulon, ulat = ds_u['lon'], ds_u['lat']
    vlon, vlat = ds_v['lon'], ds_v['lat']

    # Some models have uneven y-lengths between t/u/v — clip to common min
    if len(ulon) != len(tlon) or len(ulon) != len(vlon):
        min_len = min([len(tlon), len(ulon), len(vlon)])
        print(f'AMOC: uneven y-dim — clipping to min_len={min_len}')
        ds_t = ds_t.isel(y=slice(0, min_len))
        ds_u = ds_u.isel(y=slice(0, min_len))
        ds_v = ds_v.isel(y=slice(0, min_len))
        tlon, tlat = ds_t['lon'], ds_t['lat']
        ulon, ulat = ds_u['lon'], ds_u['lat']
        vlon, vlat = ds_v['lon'], ds_v['lat']

    grid = which_grid(tlon.values, ulon.values, vlon.values)
    print(f'AMOC: detected {grid!r}-grid')

    # Layer thickness from lev_bnds (m)
    #dz_t = ds_t['lev_bnds'].diff('bnds').isel(bnds=0)
    dz_t = ds_t.dz

    # Sigma2 from thetao, so via gsw (referenced to 2000 dbar via gsw.sigma2)
    thetao = ds_t['thetao']
    so = ds_t['so']
    lev = ds_t['lev']
    p = gsw.p_from_z(-1 * lev, tlat, geo_strf_dyn_height=0, sea_surface_geopotential=0)
    SA = gsw.SA_from_SP(so, p, tlon, tlat)
    CT = gsw.CT_from_pt(SA, thetao)
    sigma2_temp = gsw.sigma2(SA, CT)
    sigma2_temp = sigma2_temp.to_dataset(name=_AMOC_SIGNAME)[_AMOC_SIGNAME]
    sigma2_T = sigma2_temp.assign_attrs(
        {'long_name': 'Sigma referenced to 2000dbar', 'units': 'kg/m^3'}
    )

    # Bin layer thicknesses in sigma2 (gives "isopycnal thickness" — not used
    # downstream but matches the source; harmless dask graph node)
    iso_thick = histogram(
        sigma2_T, bins=[sigma_edge.values], weights=dz_t, dim=['lev'], density=False,
    )
    iso_thick = iso_thick.to_dataset(name=_AMOC_SIGNAME)[_AMOC_SIGNAME].rename(
        {_AMOC_SIGNAME + '_bin': _AMOC_SIGNAME}
    )
    iso_thick = iso_thick.assign_coords({_AMOC_SIGNAME: sigma_mid})

    # Build the horizontal volume fluxes u_e, v_e on the chosen grid type
    if use_currents:
        # Velocity path: u_e = uo * dyu * dz, v_e = vo * dxv * dz (or dxu for b-grid)
        lat_u, lon_u = ds_u['lat'], ds_u['lon']
        lat_v, lon_v = ds_v['lat'], ds_v['lon']
        if len(lat_u.shape) == 1:
            lon_u = lon_u.expand_dims(dim={'y': len(lat_u.y)})
            lat_u = lat_u.expand_dims(dim={'x': len(lon_u.x)})

        if grid in ('b', 'a'):
            htn = great_circ_dist2(lat_u, lon_u, 'x')   # x-spacing between U-centers
            hte = great_circ_dist2(lat_u, lon_u, 'y')   # y-spacing between U-centers
            dxu = (htn + htn.roll(x=-1)) / 2            # U-point centered dx
            dyu = (hte + hte.roll(y=-1)) / 2            # U-point centered dy
            u_e = ds_u['uo'] * dyu * dz_t
            v_e = ds_v['vo'] * dxu * dz_t
        elif grid == 'c':
            hte = great_circ_dist2(lat_u, lon_u, 'y')
            dyu = (hte + hte.roll(y=-1)) / 2
            htn = great_circ_dist2(lat_v, lon_v, 'x')
            dxv = (htn + htn.roll(x=-1)) / 2
            u_e = ds_u['uo'] * dyu * dz_t
            v_e = ds_v['vo'] * dxv * dz_t
        else:
            raise ValueError(f"AMOC: unsupported grid type for velocity path: {grid!r}")
    else:
        # Mass-transport path: divide by reference density to get volume flux
        ref_den = 1028
        u_e = ds_u['umo'] / ref_den
        v_e = ds_v['vmo'] / ref_den

    u_e = u_e.where(u_e < 1.0e30).fillna(0.0)
    v_e = v_e.where(v_e < 1.0e30).fillna(0.0)

    # Interpolate b- or a-grid fluxes onto c-grid (wflux_div is c-grid)
    if grid == 'b':
        u_e = 0.5 * (u_e + u_e.shift(y=1))
        v_e = 0.5 * (v_e + v_e.roll(x=1, roll_coords=False))
    elif grid == 'a':
        u_e = 0.5 * (u_e + u_e.roll(x=1, roll_coords=False))
        v_e = 0.5 * (v_e + v_e.shift(y=1))
        # Align lev to sigma2_T's lev (a-grid uses different staggering)
        u_e = u_e.interp(lev=sigma2_T.lev)
        v_e = v_e.interp(lev=sigma2_T.lev)

    # Bin u/v fluxes into sigma2 bins
    iso_uflux_temp = histogram(
        sigma2_T, bins=[sigma_edge.values], weights=u_e, dim=['lev'], density=False,
    )
    iso_uflux = iso_uflux_temp.to_dataset(name='iso_uflux').rename(
        {_AMOC_SIGNAME + '_bin': _AMOC_SIGNAME}
    ).assign_coords({_AMOC_SIGNAME: sigma_mid})

    iso_vflux_temp = histogram(
        sigma2_T, bins=[sigma_edge.values], weights=v_e, dim=['lev'], density=False,
    )
    iso_vflux = iso_vflux_temp.to_dataset(name='iso_vflux').rename(
        {_AMOC_SIGNAME + '_bin': _AMOC_SIGNAME}
    ).assign_coords({_AMOC_SIGNAME: sigma_mid})

    # Convergence to get vertical volume flux in sigma space (using c-grid
    # formulation since all fluxes have been mapped to c-grid above)
    wflux = wflux_div('c', iso_uflux['iso_uflux'], iso_vflux['iso_vflux'],
                      _AMOC_SIGNAME, sigma_edge)

    # Basin masks (Atlantic, IndoPac, Global) at first timestep
    rmaskmoc = basinmask(so.isel(time=0))

    # Cumulative meridional integration of zonally-summed wflux per basin
    moc = compute_MOC(wflux, tlat, rmaskmoc, _AMOC_LAT_MID)

    # Atlantic southern-boundary correction: find southernmost Atlantic
    # y-index, take vflux contribution south of there, integrate top-down in
    # sigma, divide by 1e6 to get Sv. Add to Atlantic MOC.
    tmp = rmaskmoc.isel(region=1).sum('x')
    atl_j = 0
    j = 0
    while atl_j == 0:
        if tmp.isel(y=j).data > 0:
            atl_j = j
        j += 1
    atl_j = atl_j - 1

    tmp = iso_vflux['iso_vflux'] * (rmaskmoc.shift(y=-1))
    tmp = tmp.chunk({'y': 1})  # deliberate chunk for high-res speed
    tmp = tmp.isel(y=atl_j, region=1).sum('x')
    moc_s = -tmp.sortby('sigma', ascending=False).cumsum('sigma').sortby('sigma', ascending=True) / 1.0e6
    moc_s['sigma'] = sigma_edge.isel({_AMOC_SIGNAME: slice(0, -1)})

    atlantic_moc = moc.isel(region=1) + moc_s
    moc = xr.concat(
        [moc.isel(region=0), atlantic_moc, moc.isel(region=2)], dim='region',
    )

    # Package output: MOC plus time_bnds from the input
    moc_ds = moc.to_dataset(name='MOC')
    moc_ds = moc_ds.assign_coords({'time': sigma2_T['time']})
   # moc_ds['time_bnds'] = ds_t['time_bnds']
    moc_ds = moc_ds.chunk(None)
    # Record which streamfunction this is, so plots label it correctly:
    # uo/vo (velocity) -> Eulerian; umo/vmo (mass transport) -> residual.
    sf_type = 'Eulerian' if use_currents else 'residual'
    moc_ds.attrs['streamfunction_type'] = sf_type
    moc_ds['MOC'].attrs['streamfunction_type'] = sf_type
    return moc_ds
