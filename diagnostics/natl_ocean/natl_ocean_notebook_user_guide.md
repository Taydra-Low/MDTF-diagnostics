# North Atlantic Ocean POD: ESNB Notebook — Setup & User Guide

Step-by-step instructions for running the `natl_ocean` POD notebook
(`diagnostics/natl_ocean/natl_ocean_esnb.ipynb`) from scratch, for someone new to the MDTF
framework. The notebook lives on the MDTF **notebook-development branch**
(`dev-notebook-transistions`), not `main`.

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

To keep the example simple and **runnable by anyone**, this notebook uses **one year** of a single
data source (`data_source = "timeslice_cmorized"`, Section 1): a CMORized copy of the CESM
(`cesm_mdtfv3`) timeslice run — CMIP-named local NetCDF (`thetao`, `so`, `uo`, `vo`, `tos`, `hfds`,
`wfo`, plus fx `areacello`/`volcello`), ~0.5 GB. It loads **through ESNB**
(`NotebookDiagnostic → CaseGroup2 → open`) like a CMIP model, and the **same files also load through
the companion intake-esm notebook** (`natl_ocean_intake_esm.ipynb`).

> **On NCAR HPC?** You can point at the local CMORized dataset and use **longer time periods** (the
> timeslice record covers 1995–2015) — see [Step 3](#step-3--get-the-data). Obs, input, and a
> prebuilt one-year bundle are also on Google Drive:
> [shared data folder](https://drive.google.com/drive/folders/1gNKwbffdD2m7d9Q7atEnbjG9pThvRU4i?usp=drive_link).

#### Why CMORized-local? (ESNB vs. intake-esm)

The two loaders differ in what storage they can open:

| | Local NetCDF | OPeNDAP (`http:`) | Cloud zarr (`gs:`) |
|---|---|---|---|
| **intake-esm** | ✅ | ✅ | ✅ (anonymous) |
| **ESNB** | ✅ | ❌ rejects `http:` | ⚠️ `gs:` path exists but omits the anon token → hangs on public data as-shipped |

So **ESNB reliably opens local files only**, while **intake-esm** can also stream OPeNDAP and cloud
data. And *native* (non-CMORized) timeslice output works with **neither** directly — both need the
`load_native_timeslice` harmonization first. *(That harmonization is specific to the **current**
CESM timeslice data, which ships as raw POP output; later iterations that provide CMORized data
upstream may not need it.)* **CMORized-local NetCDF** is the one form that loads cleanly through
**both**, which is why we use it here.

---

## Prerequisites

- A **GitHub account**.
- **SSH key** registered with GitHub (recommended). To create one and add it:
  ```bash
  ssh-keygen -t ed25519 -C "your_email@example.com"   # press enter to accept defaults
  cat ~/.ssh/id_ed25519.pub                            # copy this output
  ```

  Then paste it at GitHub → Settings → SSH and GPG keys → New SSH key. (See the MDTF git
  intro §6.10 for details: [https://mdtf-diagnostics.readthedocs.io/en/main/sphinx/dev_git_intro.html](https://mdtf-diagnostics.readthedocs.io/en/main/sphinx/dev_git_intro.html).)
- **conda** (miniconda or Anaconda) on your machine.
- Access to the input **data** — either NCAR Casper (where the data currently lives) or a
  machine you have copied the data to (see Step 3).

---

## Step 1 — Get the code (notebook-development branch)

The notebook is on the `dev-notebook-transistions` branch (note the spelling), not `main`.

**If you only want to run the notebook**, clone the main repo and check out the branch:

```bash
git clone git@github.com:NOAA-GFDL/MDTF-diagnostics.git
cd MDTF-diagnostics
git fetch origin
git checkout -b dev-notebook-transistions origin/dev-notebook-transistions
```

**If you plan to contribute changes back**, first fork the repo on GitHub (button in the
upper-right at [https://github.com/NOAA-GFDL/MDTF-diagnostics](https://github.com/NOAA-GFDL/MDTF-diagnostics)), then:

```bash
git clone git@github.com:<your_github_account>/MDTF-diagnostics.git
cd MDTF-diagnostics
git remote add upstream git@github.com:NOAA-GFDL/MDTF-diagnostics.git
git fetch upstream
git checkout -b dev-notebook-transistions upstream/dev-notebook-transistions
```

The notebook is then at:

```
diagnostics/natl_ocean/natl_ocean_esnb.ipynb
```

> Note: `dev-notebook-transistions` is the literal branch name on the remote (it contains a
> typo — "transi**sti**ons"). Use it exactly as written.

---

## Step 2 — Create the `esnb` conda environment

This POD needs more than the base esnb packages (it uses xESMF, gsw, xwmt, etc.). Create a
dedicated environment with conda for the compiled/geo packages and pip for the rest:

```bash
conda create -n esnb -c conda-forge python=3.12 xarray scipy pandas numba \
    xesmf esmpy cartopy gsw-xarray xhistogram cf_xarray matplotlib gcsfs ipykernel
conda activate esnb
pip install esnb "xwmt==0.0.3" cmip_basins intake-esm dask dask-jobqueue netcdf4 cftime
python -m ipykernel install --user --name esnb --display-name "Python (esnb)"
```

Notes:

- **`xwmt==0.0.3` is pinned on purpose.** `POD_utils.calc_wmt` calls `xwmt.swmt(...)`, a class
  that only exists in xwmt 0.0.3; xwmt ≥ 0.1.0 renamed it and changed the signature.
- The last command registers the kernel as **"Python (esnb)"**, which you select in Jupyter.

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
   (`nlat/nlon` + 2D `lon/lat`, matching CESM2's real CMIP6 files), copies the fx grid files, and
   builds the intake-esm catalog `CESM_timeslice_cmorized_001.{json,csv}`.

3. **Point the notebook at the output**: place the `.json`/`.csv` in `diagnostics/natl_ocean/` with
   the CSV `path` column set to your data.

> **A note on the fx (grid) files.** The bundle needs `areacello` (cell area) and `volcello` (cell
> volume), but the native timeslice files on Globus **don't yet include grid geometry** (no
> `TAREA`/`dz`), so those two are derived from a separate dataset — the CESM2 gx1 `Ofx` grid, whose
> geometry matches the timeslice grid. **This is already done for you:** the bundle ships the fx
> files and the notebook reads them from the catalog, so you normally don't need to do anything. If
> you *do* want them standalone, they're on ESGF — download CESM2 `areacello` / `volcello`
> (`table_id=Ofx`, `experiment_id=historical`, `variant_label=r1i1p1f1`, `grid_label=gn`) from the
> [ESGF CMIP6 search](https://esgf-node.llnl.gov/search/cmip6/) (verified available on
> `esgf-data.ucar.edu` and other nodes, plus a Globus endpoint).

### Updating catalog paths (off-NCAR)

The catalog files in `diagnostics/natl_ocean/` store **absolute** paths:

- the `path` column in `CESM_timeslice_cmorized_001.csv` (and the native
  `CESM_timeslice_ocean_monthly_001.csv` if you use `timeslice_native`)
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
| `FX_DIR`      | fx grid dir — **not needed** for `timeslice_cmorized` (the bundle carries its own fx); only used by the CESM2 sources / `timeslice_native` |
| `PBS_ACCOUNT` | your PBS account (only used if you enable the dask cluster) |

---

## Step 5 — Launch and run

1. Start Jupyter:
   - **On Casper**: use JupyterHub ([https://jupyterhub.hpc.ucar.edu](https://jupyterhub.hpc.ucar.edu)) or a `qvscode`
     interactive session; a 1-year run is small enough to run single-process (no dask).
   - **Locally**: `conda activate esnb && jupyter lab`
2. Open `diagnostics/natl_ocean/natl_ocean_esnb.ipynb` and select the **`Python (esnb)`**
   kernel (top-right kernel selector).
3. In Section 1, set `data_source` — `"timeslice_cmorized"` (default; loaded through ESNB) or
   `"timeslice_native"` (in-memory bypass) — and the date range. For a first test, use one year:
   `startdate, enddate = '1995-01-01', '1995-12-31'` (the timeslice data covers 1995–2015).
4. Run Section 1 → Section 2 (data load), then the diagnostics and plotting sections in order.

---

## Troubleshooting

Developer risks/fallbacks (xwmt pinning, WMT dim names, surface-field shapes, bias-path
fallbacks) are documented in
`diagnostics/natl_ocean/notes/2026-06-29_session_notes.md`.

Common issues:

- **Kernel "Python (esnb)" not listed** — re-run the `ipykernel install` line from Step 2.
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

**POD authors / contributors:** Elizabeth Maroon, Stephen Yeager, Taydra Low, Brendan Myers, Feng Shu

**Last updated:** 2026-06-29
