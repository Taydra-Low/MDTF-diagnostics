#!/usr/bin/env python
"""CMORize a slice of the CESM (cesm_mdtfv3) timeslice output so it loads through ESNB.

Why
---
The raw timeslice files are native POP output (TEMP/SALT/UVEL on z_t in cm, TLONG/TLAT,
no tos, native units). ESNB opens local files with xr.open_mfdataset and expects the
CMIP variable name in-file, so native files can't be used directly. POD_utils.
load_native_timeslice() already performs the full native->CMIP harmonization in-memory
(rename, unit fixes, time recenter, the NCAR/POP2 frazil recipe for hfds/wfo, deriving
tos). This script runs that harmonizer over a chosen window, applies light CF finishing,
writes the result as CMIP-named NetCDF, copies the matching fx files, and builds an
intake-esm catalog. The written bundle then loads through ESNB like a CMORized model.

Grid convention: native POP already uses nlat/nlon dims with 2D lon/lat coords, which
matches CESM2's real CMIP6 files (CESM2 ships nlat/nlon, a deviation from the canonical
CMOR i/j; we intentionally match CESM2 so the POD treats timeslice identically). Default
is one year (1995). The MDTF framework's own CMOR tooling is not used (it does not handle
the 2D surface fields here); this is a self-contained xarray writer.

Usage
-----
    python cmorize_timeslice.py --dry-run                 # thetao+tos, 2 months, temp dir
    python cmorize_timeslice.py --year 1995               # full one-year bundle (~1.5 GB)
    python cmorize_timeslice.py --start 1995-01-01 --end 1995-12-31 --outdir /path
"""
import os
import sys
import json
import shutil
import argparse
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import POD_utils

# CMIP standard_name / long_name for the vars load_native_timeslice returns.
# Only ADD these when absent; units are left as the harmonizer set them.
CMIP_META = {
    "thetao": ("sea_water_potential_temperature",         "Sea Water Potential Temperature"),
    "so":     ("sea_water_salinity",                      "Sea Water Salinity"),
    "uo":     ("sea_water_x_velocity",                    "Sea Water X Velocity"),
    "vo":     ("sea_water_y_velocity",                    "Sea Water Y Velocity"),
    "hfds":   ("surface_downward_heat_flux_in_sea_water", "Downward Heat Flux at Sea Water Surface"),
    "wfo":    ("water_flux_into_sea_water",               "Water Flux into Sea Water"),
    "tos":    ("sea_surface_temperature",                 "Sea Surface Temperature"),
}
WRITE_ORDER = ("thetao", "so", "uo", "vo", "hfds", "wfo", "tos")

# fx files (CESM2 gx1 Ofx; the timeslice gx1 grid matches). Copied into the bundle
# so it is self-contained. Relative to FX_DIR.
FX_FILES = {
    "areacello": "areacello/gn/v20190308/areacello_Ofx_CESM2_historical_r1i1p1f1_gn.nc",
    "volcello":  "volcello/gn/v20190308/volcello_Ofx_CESM2_historical_r1i1p1f1_gn.nc",
}
FX_DIR_DEFAULT = "/glade/collections/cmip/CMIP6/CMIP/NCAR/CESM2/historical/r1i1p1f1/Ofx"


def cmor_finish(ds, var):
    """Light CF finishing: lev cm->m + CF attrs on coords and the variable."""
    ds = ds.copy()
    if "lev" in ds.coords:
        lev = ds["lev"]
        if float(lev.max()) > 8000:                 # POP z_t still in cm
            lev = lev / 100.0
        lev.attrs.update({"units": "m", "standard_name": "depth",
                          "positive": "down", "axis": "Z", "long_name": "depth"})
        ds = ds.assign_coords(lev=lev)
    if "lon" in ds.coords:
        ds["lon"].attrs.setdefault("standard_name", "longitude")
        ds["lon"].attrs.setdefault("units", "degrees_east")
    if "lat" in ds.coords:
        ds["lat"].attrs.setdefault("standard_name", "latitude")
        ds["lat"].attrs.setdefault("units", "degrees_north")
    if "time" in ds.coords:
        ds["time"].attrs.setdefault("standard_name", "time")
        ds["time"].attrs.setdefault("axis", "T")
    if var in ds.data_vars and var in CMIP_META:
        sn, ln = CMIP_META[var]
        ds[var].attrs.setdefault("standard_name", sn)
        ds[var].attrs.setdefault("long_name", ln)
    return ds


def time_range_str(ds):
    t = ds["time"]
    def fmt(x):
        try:
            return f"{x.year:04d}{x.month:02d}"
        except AttributeError:
            ts = pd.Timestamp(x); return f"{ts.year:04d}{ts.month:02d}"
    return f"{fmt(t.values[0])}-{fmt(t.values[-1])}"


def out_filename(var, ds):
    return f"{var}_Omon_CESM-mdtfv3-timeslice_timeslice_r1i1p1f1_gn_{time_range_str(ds)}.nc"


def write_var(ds, var, outdir, compress=True):
    ds = cmor_finish(ds, var)
    ds.attrs.update({
        "Conventions": "CF-1.7", "source_id": "CESM-mdtfv3-timeslice",
        "experiment_id": "timeslice", "variable_id": var, "table_id": "Omon",
        "grid_label": "gn", "frequency": "mon", "realm": "ocean",
        "member_id": "r1i1p1f1",
        "cmorized_by": "cmorize_timeslice.py (custom xarray writer)",
    })
    path = Path(outdir) / out_filename(var, ds)
    enc = {var: {"zlib": True, "complevel": 1, "_FillValue": 1.0e20}} if compress else {}
    ds.to_netcdf(path, encoding=enc)
    return str(path), time_range_str(ds)


def build_catalog(rows, native_csv, outdir, base="CESM_timeslice_cmorized_001"):
    """rows: list of dicts (variable_id, path, table_id, time_range)."""
    src = pd.read_csv(native_csv)
    tmpl = src.iloc[0]
    out_rows = []
    for r in rows:
        srow = src[src["variable_id"] == r["variable_id"]]
        cm = (srow["cell_methods"].iloc[0] if not srow.empty
              else tmpl.get("cell_methods", "area: mean time: mean"))
        out_rows.append({
            "activity_id":   tmpl.get("activity_id", "CESM"),
            "institution_id": tmpl.get("institution_id", "NCAR"),
            "source_id":     "CESM-mdtfv3-timeslice",
            "experiment_id": "timeslice", "frequency": "mon", "realm": "ocean",
            "variable_id":   r["variable_id"], "native_variable_name": "",
            "member_id":     tmpl.get("member_id", "r1i1p1f1"),
            "table_id":      r["table_id"], "grid_label": "gn",
            "cell_methods":  cm, "chunk_freq": tmpl.get("chunk_freq", ""),
            "time_range":    r["time_range"], "path": r["path"],
            "version":       tmpl.get("version", "v1"),
        })
    cat = pd.DataFrame(out_rows, columns=list(src.columns))
    csv_path = Path(outdir) / f"{base}.csv"
    cat.to_csv(csv_path, index=False)
    with open(Path(native_csv).with_suffix(".json")) as f:
        meta = json.load(f)
    meta["id"] = base
    meta["catalog_file"] = str(csv_path)
    meta["description"] = "CMORized CESM mdtfv3 timeslice ocean monthly (custom writer)"
    json_path = Path(outdir) / f"{base}.json"
    with open(json_path, "w") as f:
        json.dump(meta, f, indent=2)
    return str(csv_path), str(json_path)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalog", default="CESM_timeslice_ocean_monthly_001.json")
    ap.add_argument("--outdir", default="/glade/derecho/scratch/taydral/cesm_mdtfv3_timeslice_cmorized_1yr")
    ap.add_argument("--year", type=int, default=None, help="convenience: sets --start/--end to that calendar year")
    ap.add_argument("--start", default="1995-01-01")
    ap.add_argument("--end",   default="1995-12-31")
    ap.add_argument("--fx-dir", default=FX_DIR_DEFAULT)
    ap.add_argument("--no-fx", action="store_true", help="do not copy fx (areacello/volcello) into the bundle")
    ap.add_argument("--no-compress", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    if args.year is not None:
        args.start, args.end = f"{args.year}-01-01", f"{args.year}-12-31"

    catalog_json = str(Path(args.catalog).resolve()) if os.path.exists(args.catalog) \
        else str(Path(__file__).resolve().parent / args.catalog)
    native_csv = str(Path(catalog_json).with_suffix(".csv"))
    print(f"native catalog: {catalog_json}")
    print(f"window: {args.start} .. {args.end}")

    datasets = POD_utils.load_native_timeslice(
        catalog_json, varnames=None, startdate=args.start, enddate=args.end)
    by_var = {list(d.data_vars)[0]: d for d in datasets}
    print(f"harmonized variables: {sorted(by_var)}")
    print(f"months in window: {by_var['thetao'].sizes.get('time')}")

    if args.dry_run:
        tmp = tempfile.mkdtemp(prefix="cmor_dryrun_")
        ok = True
        for var in ("thetao", "tos"):
            sub = by_var[var].isel(time=slice(0, 2))
            p, _ = write_var(sub, var, tmp, compress=not args.no_compress)
            re = xr.open_dataset(p)
            has = var in re.data_vars
            lev_ok = ("lev" not in re.coords) or (float(re["lev"].max()) < 8000)
            grid_ok = ("nlat" in re.dims and "nlon" in re.dims)
            print(f"[dry-run] {var}: {os.path.basename(p)} {os.path.getsize(p)/1e6:.1f} MB "
                  f"| in-file var={has} | lev_in_m={lev_ok} | nlat/nlon={grid_ok} | dims={dict(re.sizes)}")
            ok = ok and has and lev_ok and grid_ok
            re.close()
        shutil.rmtree(tmp, ignore_errors=True)
        print("[dry-run] PASS" if ok else "[dry-run] FAIL")
        sys.exit(0 if ok else 1)

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    print(f"outdir: {outdir}")
    rows = []
    for var in WRITE_ORDER:
        if var not in by_var:
            print(f"  SKIP {var} (not produced)"); continue
        print(f"  writing {var} ...", flush=True)
        p, tr = write_var(by_var[var], var, outdir, compress=not args.no_compress)
        rows.append({"variable_id": var, "path": p, "table_id": "Omon", "time_range": tr})
        print(f"    -> {os.path.basename(p)} ({os.path.getsize(p)/1e9:.2f} GB)")

    if not args.no_fx:
        fxbase = Path(args.fx_dir)
        for var, rel in FX_FILES.items():
            src = fxbase / rel
            dst = outdir / f"{var}_Ofx_CESM-mdtfv3-timeslice_timeslice_r1i1p1f1_gn.nc"
            if src.exists():
                shutil.copyfile(src, dst)
                rows.append({"variable_id": var, "path": str(dst), "table_id": "Ofx", "time_range": ""})
                print(f"    fx -> {dst.name} ({os.path.getsize(dst)/1e6:.1f} MB)")
            else:
                print(f"    fx MISSING: {src} (skipped)")

    csv_path, json_path = build_catalog(rows, native_csv, outdir)
    total = sum(os.path.getsize(r["path"]) for r in rows) / 1e9
    print(f"catalog CSV : {csv_path}")
    print(f"catalog JSON: {json_path}")
    print(f"done: {len(rows)} files, {total:.2f} GB total")


if __name__ == "__main__":
    main()
