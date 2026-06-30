# natl_ocean POD — session notes (2026-06-29)

Work on `example_notebooks/mdtf.natl_ocean.esnb.ipynb` and the `natl_ocean` POD: wiring up
the CESM timeslice data source, de-hardcoding the notebook for use as a template, filling
out the docs, and capturing developer risks/fallbacks (moved here out of the notebook).

## Summary of changes

### 1. Loading errors fixed (CESM timeslice)
The second data source was mislabeled "GFDL_timeslice"; it is actually **CESM** timeslice
output (`cesm_mdtfv3_timeslice`, native POP on the gx1 grid). Three successive `.open()`
failures were diagnosed and resolved:
- **`OSError: no files to open`** — the requested date range was outside the catalog's
  coverage (`19950101–20150101`). Default dates set inside that window.
- **`FileNotFoundError` on the catalog CSV** — the GFDL catalog JSON stored a bare relative
  `catalog_file`; CESM2's used an absolute `file://` path. Made the notebook resolve a
  relative `catalog_file` against the catalog JSON's directory, and set the catalog JSON to
  an absolute path.
- **`KeyError: "No variable named 'thetao'"`** — ESNB opens by CMIP `variable_id` and
  requires that name in-file, but the native POP files carry native names (`TEMP`, `SALT`,
  …). Root cause: the timeslice data is **not CMORized**. Fix: harmonize in-memory (below).

### 2. Rename GFDL → CESM
- Catalog files `GFDL_timeslice_ocean_monthly_001.{json,csv}` →
  `CESM_timeslice_ocean_monthly_001.{json,csv}`; column `gfdl_variable_name` →
  `native_variable_name`; `institution_id` GFDL → NCAR; removed stale `.csv.bak`.
- Notebook: `data_source = "CESM_timeslice"`, all branches, comments, and markdown updated.
- Left the title *"Unified MDTF/GFDL/NCAR Analysis Notebook Template"* unchanged — "GFDL/NCAR"
  there refers to the two MDTF host centers, not this dataset.

### 3. In-notebook harmonization (`load_native_timeslice`)
Section 2 now branches on `data_source`: CESM2 loads via ESNB; **CESM_timeslice bypasses
ESNB** and harmonizes native POP files to CMIP in-memory. Driven by a `mapping` dict
`{cmip: native}` built from the catalog's `native_variable_name` column (MCS-notebook style).
Steps:
- rename native variable → CMIP id; coords `TLONG/TLAT`, `ULONG/ULAT` → `lon/lat`, `z_t` → `lev`;
- normalize units (`gram/kilogram`→`g/kg`, `watt/m^2`/`Watts/meter^2`→`W m-2`,
  `kg/m^2/s`→`kg m-2 s-1`, velocities `cm/s`→`m/s` ×0.01);
- recenter POP end-of-interval monthly time stamps to the interval midpoint (via `time_bound`);
- derive `tos` from surface `thetao` (native output has no SST);
- open lazily (`chunks={'time': 1}`) so unit math doesn't materialize full arrays.
Both branches produce `loaded_datasets`, consumed by the existing per-variable dispatch cell.

### 4. CESM/POP2 frazil surface-flux recipe (CESM-only)
Corrected the `hfds` mapping from `QFLUX` (frazil ice-formation heat flux) to **`SHF`** (total
surface heat flux). The catalog now maps `SHF→hfds`, `QFLUX→hfsifrazil`, `SFWF→wfo`. The
loader applies the standard NCAR recipe (mirrors `sub2sub_HR/wmt.py::_op_frazil_correction`
and the OSNAP NCAR1/NCAR10 path):
- `hfds = SHF + QFLUX`
- `wfo  = SFWF − QFLUX / (latent_heat_fusion / 1e4)`, with
  `latent_heat_fusion = 3.337e9 erg/g` (read from CESM FOSI POP output → 3.337e5 J/kg).
This is applied **only to native CESM/POP output** — CMORized models (CESM2 via ESNB) skip it.

### 5. De-hardcoding for portability (template)
Single top-of-notebook config block via `os.environ.setdefault`: `CODE_ROOT`, `WORK_DIR`,
`OBS_DATA`, `FX_DIR`, `PBS_ACCOUNT`, `CONDA_ROOT` (derived from `CONDA_EXE`). All downstream
paths reference these:
- fixed the `OBS_DATA` bug (was commented out but used by the obs-WMT cells);
- conda root / env root derived rather than literal;
- fx dir via `FX_DIR`; obs file and OMIP benchmark built from `$OBS_DATA`
  (`omip2.cycle1.1989_2018.0-200m.mld_sic_t_s_sigma.nc` now lives in the obs_data dir);
- PBS account via `PBS_ACCOUNT` (dask cell parameterized, not run);
- removed resolved author TODO notes.

### 6. Documentation (`doc/natl_ocean.rst`)
Filled out from the MDTF skeleton: synopsis, functionality (added WMT + AMOC paths),
full required-variable table, note on the CESM-only frazil recipe, references (kept
`ref-Maloney`; science papers listed as *in preparation* with a TODO for the DOI — papers
are still in progress). Removed a documented function that does not exist (`reorg_by_cycle`).

## Verification (1 year, 1995; single-process, no dask)
- Loader: 7 catalog vars + derived `tos`; `hfsifrazil` consumed; units harmonized; lev→m;
  12 months selected (no POP off-by-one); `tos == surface thetao`.
- Frazil recipe matched a manual computation **exactly** at both a zero-frazil point and a
  high-frazil point (QFLUX = 13.2 W/m²): `hfds = SHF + QFLUX`, `wfo = SFWF − QFLUX/latfus`.
- Full path loader → dispatch → bias merge → `preprocess_coords` → σ₀ runs clean (σ₀ finite).
- CESM2 ESNB resolve still works (no regression); all notebook cells compile; notebook valid.

## Open items / things to confirm
- **Flux sign conventions** (CMIP positive-down) not independently verified against POP —
  check if WMT magnitudes look off.
- Freshwater forcing for the timeslice now uses `wfo` (frazil-corrected SFWF); the dispatch
  prefers `wfo` over `vsf`.
- `doc/natl_ocean.rst` references: add DOIs once the science papers are published.

---

## Developer risks / fallbacks (moved out of the notebook)

**Bias-plot path (Section 4 first half):**
- **`POD_utils.preprocess_coords`** may rename `nlat`/`nlon` to other names. If downstream
  code errors after Step 1, fall back to `ds_target = ds_bias.copy()` (skip the preprocess call).
- **`compute_zavg`** requires `dz` to have `lev` in meters and broadcast against the variable.
  If `*_zavg` arrays come back with NaNs or wrong shape, recheck the cm→m conversion and
  `dz.assign_coords(lev=...)`.
- **`SpatialPlot_climo_bias`** expects `ds_model` to carry variables named `thetao_zavg`,
  `so_zavg`, `sigma0_zavg`, `mld`. The OMIP file's actual var names may differ — inspect
  `ds_model.data_vars` after loading and rename with `ds_model = ds_model.rename({...})` if needed.
- **`ScatterPlot_Error`** assumes the bias dataset returned by `SpatialPlot_climo_bias` contains
  the `*_bias` fields. If it errors, inspect `ds_t200.data_vars` after the first plot and adjust.

**WMT path (Section 3 Steps 6–7 + Section 4 second half):**
- **xwmt version pin**: `POD_utils.calc_wmt` calls `xwmt.swmt(...)` — that class only exists
  in xwmt 0.0.3. xwmt ≥ 0.1.0 replaced it with `WaterMassTransformations`, which has a different
  signature. The esnb env is pinned to xwmt 0.0.3 for this reason.
- **Surface-fields shape**: Step 6 builds `ds_target_wmt` independently of the bias-path
  `ds_target` (which has been regridded + climatology'd). WMT needs the native-grid time series,
  not the climatology.
- **`wmt_plot_maps` dim names**: the third arg expects `[time_coord, lon_coord, lat_coord]`
  matching the result of `compute_wmt(..., 'MAPS')`. We pass `['time', 'x', 'y']` since `wmt_preproc`
  renames spatial dims. If it errors, inspect `ds_wmt_maps.dims` and adjust the list.
- **xwmt `vsf` handling**: `wmt_preproc` converts `vsf` to `wfo` (the CMIP wfo would have the
  opposite sign convention). We feed `vsf` so this branch is exercised.
