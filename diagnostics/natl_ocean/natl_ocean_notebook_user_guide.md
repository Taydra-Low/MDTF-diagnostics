# North Atlantic Ocean POD: ESNB Notebook — Setup & User Guide

Step-by-step instructions for running the `natl_ocean` POD notebook
(`diagnostics/natl_ocean/natl_ocean_esnb.ipynb`) from scratch, for someone new to the MDTF
framework. The notebook lives on the **`feature/taydra-dev`** branch of the
**Taydra-Low/MDTF-diagnostics** fork ([github.com/Taydra-Low/MDTF-diagnostics](https://github.com/Taydra-Low/MDTF-diagnostics)).

> **Developing with Claude?** For how this POD was built with Claude Code on NCAR HPC — the CLAUDE.md
> setup, HPC/compute-node conventions, anti-hallucination rules, and a reusable prompt library — see
> the companion [`claude_workflow_guide.md`](./claude_workflow_guide.md).

---

## Overview

The North Atlantic Ocean POD (`natl_ocean`) evaluates how well a model represents the
subtropical-to-subpolar North Atlantic ocean state, compared against
observation-based data. It computes:

- Sea-surface and upper-200 m temperature/salinity climatologies and biases
- Mixed layer depth (MLD) and potential density (σ₀)
- Surface heat and freshwater flux fields
- Surface-forced water-mass transformation (WMT) in σ₂ coordinates, by region
- Atlantic Meridional Overturning Circulation (AMOC) streamfunction in σ₂ coordinates
- A WMT + AMOC-at-45°N compensation diagnostic

### Data source: CMORized CESM timeslice

To keep the example simple and **runnable by anyone**, this notebook uses **one year** of a single,
fixed data source: a CMORized copy of the CESM (`cesm_mdtfv3`) timeslice run — CMIP-named local
NetCDF (`thetao`, `so`, `uo`, `vo`, `tos`, `hfds`, `wfo`/`vsf`, plus fx `areacello`/`volcello`),
~0.5 GB. There is no data-source switch: Section 1 points at the bundle's catalog
(`CESM_timeslice_cmorized_001.json`); the editable knobs are `mode`, `machine`, and `dask` (see
[Step 4](#step-4--configure-paths)). It loads **through ESNB**
(`NotebookDiagnostic → CaseGroup2 → open`) like a CMIP model, and the **same files also load through
the companion intake-esm notebook** (`natl_ocean_intake_esm.ipynb`). The Intake-ESM notebook is still a work in progress.

> **On NCAR HPC?** You can point at the local CMORized dataset and use **longer time periods** (the
> timeslice record covers 1995–2015) — see [Step 3](#step-3--get-the-data). Obs, input, and a
> prebuilt one-year bundle are also on Google Drive:
> [shared data folder](https://drive.google.com/drive/folders/1gNKwbffdD2m7d9Q7atEnbjG9pThvRU4i?usp=drive_link).

#### Why CMORized-local? (ESNB vs. Intake-ESM)

The two loaders differ in what storage they can open:

|                      | Local NetCDF | OPeNDAP (`http:`) | Cloud zarr (`gs:`)                                                                |
| -------------------- | ------------ | ------------------- | ----------------------------------------------------------------------------------- |
| **intake-esm** | ✅           | ✅                  | ✅ (anonymous)                                                                      |
| **ESNB**       | ✅           | ❌ rejects`http:` | ⚠️`gs:` path exists but omits the anon token → hangs on public data as-shipped |

So **ESNB reliably opens local files only**, while **intake-esm** can also stream OPeNDAP and cloud
data. And *native* (non-CMORized) timeslice output works with **neither** directly — both need the
`load_native_timeslice` harmonization first. *(That harmonization is specific to the **current**
CESM timeslice data, which ships as raw POP output; later iterations that provide CMORized data
upstream may not need it.)* **CMORized-local NetCDF** is the one form that loads cleanly through
**both**, which is why we use it here.

---

## Prerequisites

- **git** (to clone the repo over HTTPS — no GitHub account or SSH key required just to run the
  example).
- **conda** (miniconda or Anaconda) on your machine.
- Access to the input **data** — either NCAR Casper (where the data currently lives) or a
  machine you have copied the data to (see Step 3).

---

## Step 1 — Get the code

The POD lives on the `feature/taydra-dev` branch of the Taydra-Low fork. Clone it over HTTPS and
check out the branch:

```bash
git clone https://github.com/Taydra-Low/MDTF-diagnostics.git
cd MDTF-diagnostics
git checkout feature/taydra-dev
```

The notebook is then at:

```
diagnostics/natl_ocean/natl_ocean_esnb.ipynb
```

---

## Step 2 — Create the `mdtf_nb` conda environment

The notebook runs in a dedicated `mdtf_nb` environment (its Jupyter kernel is named `mdtf_nb`). It
needs more than the base esnb packages (it uses xESMF, gsw, xwmt, regionmask, etc.). Create it with
conda for the compiled/geo packages and pip for the rest:

```bash
conda create -n mdtf_nb -c conda-forge python=3.12 xarray scipy pandas numba \
    xesmf esmpy cartopy gsw-xarray xhistogram cf_xarray regionmask matplotlib gcsfs ipykernel
conda activate mdtf_nb
pip install esnb "xwmt==0.0.3" cmip_basins intake-esm dask dask-jobqueue netcdf4 cftime
python -m ipykernel install --user --name mdtf_nb --display-name "mdtf_nb"
```

Notes:

- **`xwmt==0.0.3` is pinned on purpose.** `POD_utils.calc_wmt` calls `xwmt.swmt(...)`, a class
  that only exists in xwmt 0.0.3; xwmt ≥ 0.1.0 renamed it and changed the signature.
- The last command registers the kernel as **`mdtf_nb`**, which you select in Jupyter.

---

## Step 3 — Get the data

The notebook needs two sets of files: **reference/obs data** and the **CMORized timeslice model
data** (which bundles its own fx grid files). On NCAR Casper they already exist at the paths below;
off-NCAR, get them from Google Drive or rebuild the model data yourself (below).

**a) Reference / observation data** (point `OBS_DATA` at this directory):

- `obs_1x1.nc`
- `omip2.cycle1.1989_2018.0-200m.mld_sic_t_s_sigma.nc`
- `obs.maps_freq.sigma2.1982-2009_decomp_mean_1x1.nc`
- `obs_wmt_sigma2_1982-2009_spna_decomp_mean.nc`

On Casper: `/glade/work/taydral/MDTF/inputdata/obs_data/natl_ocean/`. Off-NCAR: download from the
[shared data folder (Google Drive)](https://drive.google.com/drive/folders/1gNKwbffdD2m7d9Q7atEnbjG9pThvRU4i?usp=drive_link).

**b) Model data — the CMORized CESM timeslice**

The notebook needs the CMORized files (`CESM_timeslice_cmorized_001.{json,csv}` + the `*.nc`,
including `areacello`/`volcello`).

> **Note:** these notebooks and the shared Globus/Drive data are a **work in progress** — this
> CMORize step may not be needed in the future as the workflow matures.

**Option A — download the prebuilt one-year bundle (easiest; recommended if you're new to Globus).**
~0.5 GB, from Google Drive:
[shared data folder](https://drive.google.com/drive/folders/1gNKwbffdD2m7d9Q7atEnbjG9pThvRU4i?usp=drive_link).
It already includes the fx grid files, so **nothing else is required**. Put the `.json`/`.csv` in
`diagnostics/natl_ocean/` and set the CSV `path` column to your local copy. On Casper the bundle is
at `/glade/work/taydral/MDTF/inputdata/cesm_mdtfv3_timeslice_cmorized_1yr/`.

**Option B — build it yourself from the native data (more data-intensive; lets you use a larger
time range).** Do this if you want more than one year.

Note: This option is untested on other machines.

1. **Get the native timeslice files via Globus.** These native POP files are the CMORization input,
   including **`QFLUX`** — the frazil ice-formation heat flux the surface-flux recipe needs (the
   latent-heat-of-fusion constant is already in the code, so no extra conversion file is required):
   `cesm_mdtfv3_timeslice.{TEMP,SALT,SHF,QFLUX,SFWF,UVEL,VVEL}.mon.nc`.
   Globus endpoint:
   [native timeslice data](https://app.globus.org/file-manager?origin_id=200c3a02-0c49-4e3c-ad24-4a24db9b1c2d&origin_path=%2F).
   On Casper: `/glade/campaign/cgd/amp/bundy/mdtf/cesm_mdtfv3_timeslice_public/ocn/mon/`. Point the
   `path` column of `diagnostics/natl_ocean/CESM_timeslice_ocean_monthly_001.csv` at them.
   *(Confirm redistribution terms with the data owner before re-sharing native files.)*
2. **Run the CMORizer** (`diagnostics/natl_ocean/cmorize_timeslice.py`, written mainly for **CESM**
   POP output):

   ```bash
   cd diagnostics/natl_ocean
   python cmorize_timeslice.py --dry-run                       # quick check: 2 months, temp dir
   python cmorize_timeslice.py --year 1995 --outdir <OUT>      # one year, ~0.5 GB
   # full record (1995–2015): ~29 GB native in -> ~13 GB CMORized (compressed) out
   python cmorize_timeslice.py --start 1995-01-01 --end 2014-12-31 --outdir <OUT>
   ```

   Run multi-year builds **under PBS**, not on the login node. The script reuses the same
   `load_native_timeslice` harmonization (CMIP renames, unit fixes, depth cm→m, derive `tos`,
   CESM/POP2 frazil recipe `hfds = SHF + QFLUX`, `wfo = SFWF − QFLUX/L`), writes CMIP-named NetCDF
   (`nlat/nlon` + 2D `lon/lat`, matching CESM2's real CMIP6 files), and builds the intake-esm catalog
   `CESM_timeslice_cmorized_001.{json,csv}`. **fx (`areacello`/`volcello`) are not written by
   default** — see the fx note below.
3. **Point the notebook at the output**: place the `.json`/`.csv` in `diagnostics/natl_ocean/` with
   the CSV `path` column set to your data.
4. **Add the fx (grid) files** (see the fx note below): download `areacello`/`volcello` from the
   shared Google Drive folder and add two rows to the catalog CSV
   (`variable_id`, `table_id=Ofx`, `path=<your file>`).

> **A note on the fx (grid) files.** The notebook needs `areacello` (cell area) and `volcello` (cell
> volume), but the native timeslice files on Globus **don't include grid geometry** (no `TAREA`/`dz`),
> so those two come from a separate dataset — the CESM2 gx1 `Ofx` grid, whose geometry matches the
> timeslice grid. The CMORizer **does not bundle them by default**: download `areacello`/`volcello`
> from the [shared Google Drive folder](https://drive.google.com/drive/folders/1gNKwbffdD2m7d9Q7atEnbjG9pThvRU4i?usp=drive_link)
> and add a catalog row for each (`table_id=Ofx`, `path` = your local copy). *On NCAR only*, you can
> instead let the CMORizer copy them from the glade CESM2 `Ofx` grid with
> `python cmorize_timeslice.py ... --copy-fx`. They're also on ESGF — CESM2 `areacello`/`volcello`
> (`table_id=Ofx`, `experiment_id=historical`, `variant_label=r1i1p1f1`, `grid_label=gn`) via the
> [ESGF CMIP6 search](https://esgf-node.llnl.gov/search/cmip6/) (verified on `esgf-data.ucar.edu`
> and other nodes, plus a Globus endpoint).

### Updating catalog paths (off-NCAR)

The catalog files in `diagnostics/natl_ocean/` store **absolute** paths:

- the `path` column in `CESM_timeslice_cmorized_001.csv`
- the `catalog_file` field in the matching `.json`

After copying the model data locally, edit those `path` entries to point at your local copies
(or regenerate the catalog with MDTF's `catalog_builder` tool). The notebook resolves a
relative `catalog_file` against the catalog JSON's own directory, so a relative path works too.

---

## Step 4 — Configure paths

All site-specific paths live in one **USER CONFIGURATION** cell near the top of the notebook
(Section 1). Edit the `CHANGE ME` lines, or export the same variables in your shell before
launching Jupyter (the cell uses `setdefault`, so exported values win):

| Variable        | What to set it to                                           |
| --------------- | ----------------------------------------------------------- |
| `CODE_ROOT`   | your`MDTF-diagnostics` repo directory                     |
| `WORK_DIR`    | output / work directory (created if missing)                |
| `OBS_DATA`    | the reference/obs directory from Step 3a                    |
| `PBS_ACCOUNT` | your PBS account (only used if you enable the dask cluster) |

(`CONDA_ROOT` is auto-derived from the active env — you don't set it. The notebook does **not** use an
`FX_DIR`: the CMORized bundle carries its own fx grid via the catalog. `FX_DIR` only matters if you
build the bundle yourself with `cmorize_timeslice.py` — see Step 3, Option B.)

A second **USER SETTINGS** cell (just below) holds run-mode toggles, normally left at their defaults:

| Setting     | Default           | Change it when…                                                             |
| ----------- | ----------------- | ---------------------------------------------------------------------------- |
| `mode`    | `"interactive"` | leave as-is;`"prod"` reads `settings.jsonc` when the framework drives it |
| `machine` | `"casper"`      | set to`"local"` if you're not on Casper                                    |
| `dask`    | `False`         | `True` to spin up a Casper PBS dask cluster (not needed for a 1-year run)  |
| `verbose` | `True`          | `False` to quiet the load/inventory prints                                 |

---

## Step 5 — Launch and run

1. Start Jupyter:
   - **On Casper**: use JupyterHub ([https://jupyterhub.hpc.ucar.edu](https://jupyterhub.hpc.ucar.edu)) or a `qvscode`
     interactive session; a 1-year run is small enough to run single-process (no dask).
   - **Locally**: `conda activate mdtf_nb && jupyter lab`
2. Open `diagnostics/natl_ocean/natl_ocean_esnb.ipynb` and select the **`mdtf_nb`**
   kernel (top-right kernel selector).
3. Section 1 already points at the CMORized timeslice catalog (loaded through ESNB) — there's no
   data-source to choose. Just set the date range; for a first test, use one year:
   `startdate, enddate = '1995-01-01', '1995-12-31'` (the timeslice data covers 1995–2015, but the data provided is from the start and end date provided above).
4. Run Section 1 → Section 2 (data load), then the diagnostics and plotting sections in order.

---

## Troubleshooting

Developer risks/fallbacks (xwmt pinning, WMT dim names, surface-field shapes, bias-path fallbacks)
are documented inline in `diagnostics/natl_ocean/POD_utils.py`.

Common issues:

- **Kernel `mdtf_nb` not listed** — re-run the `ipykernel install` line from Step 2.
- **`ModuleNotFoundError` / `xwmt.swmt` missing** — confirm `xwmt==0.0.3` (`pip show xwmt`).
- **`OSError: no files to open`** — the requested date range is outside the catalog's coverage
  (timeslice = 1995–2015), or the catalog `path` entries don't point at real files. Check the
  date range and the catalog CSV paths.
- **`FileNotFoundError` on the catalog CSV** — the catalog JSON's `catalog_file` points
  somewhere that doesn't exist; use an absolute `file://` path or a relative filename next to
  the JSON.
- **Coordinate mismatch after Step 1** — `POD_utils.preprocess_coords` standardizes coord
  names; if a downstream cell errors, inspect `ds_target.dims` / `.data_vars`.

---

## References

- MDTF-diagnostics docs: [https://mdtf-diagnostics.readthedocs.io/](https://mdtf-diagnostics.readthedocs.io/)
- Git setup guide: [https://mdtf-diagnostics.readthedocs.io/en/main/sphinx/dev_git_intro.html](https://mdtf-diagnostics.readthedocs.io/en/main/sphinx/dev_git_intro.html)
- POD detailed description: `diagnostics/natl_ocean/doc/natl_ocean.rst`
- ESNB: [https://github.com/jkrasting/esnb](https://github.com/jkrasting/esnb)
- intake-esm: [https://intake-esm.readthedocs.io/](https://intake-esm.readthedocs.io/)

---

**POD authors / contributors:** Taydra Low, Elizabeth Maroon, Stephen Yeager, Brendan Myers,
Teagan King, Feng Zhu

**Last updated:** 2026-07-06
