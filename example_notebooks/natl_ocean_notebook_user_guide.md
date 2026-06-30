# North Atlantic Ocean POD: ESNB Notebook — Setup & User Guide

Step-by-step instructions for running the `natl_ocean` POD notebook
(`example_notebooks/mdtf.natl_ocean.esnb.ipynb`) from scratch, for someone new to the MDTF
framework. The notebook lives on the MDTF **notebook-development branch**
(`dev-notebook-transistions`), not `main`.

---

## Overview

The North Atlantic Ocean POD (`natl_ocean`) evaluates how well a model represents the
subtropical-to-subpolar North Atlantic ocean state, compared against a 1°×1°
observation-based reference. It computes:

- Sea-surface and upper-200 m temperature/salinity climatologies and biases
- Mixed layer depth (MLD) and potential density (σ₀)
- Surface heat and freshwater flux fields
- Surface-forced water-mass transformation (WMT) in σ₂ coordinates, by region
- Atlantic Meridional Overturning Circulation (AMOC) streamfunction in σ₂ coordinates
- A WMT + AMOC-at-45°N compensation diagnostic

It supports two data sources, selected by the `data_source` variable in Section 1:

- **`CESM2`** — CMIP-conformant (CMORized) output, loaded through the ESNB framework.
- **`CESM_timeslice`** — native CESM/POP output, which is *not* CMORized; the notebook
  harmonizes it to CMIP variables/coordinates/units in-memory and applies the CESM/POP2
  frazil surface-flux recipe.

---

## Prerequisites

- A **GitHub account**.
- **SSH key** registered with GitHub (recommended). To create one and add it:
  ```bash
  ssh-keygen -t ed25519 -C "your_email@example.com"   # press enter to accept defaults
  cat ~/.ssh/id_ed25519.pub                            # copy this output
  ```
  Then paste it at GitHub → Settings → SSH and GPG keys → New SSH key. (See the MDTF git
  intro §6.10 for details: <https://mdtf-diagnostics.readthedocs.io/en/main/sphinx/dev_git_intro.html>.)
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
upper-right at <https://github.com/NOAA-GFDL/MDTF-diagnostics>), then:

```bash
git clone git@github.com:<your_github_account>/MDTF-diagnostics.git
cd MDTF-diagnostics
git remote add upstream git@github.com:NOAA-GFDL/MDTF-diagnostics.git
git fetch upstream
git checkout -b dev-notebook-transistions upstream/dev-notebook-transistions
```

The notebook is then at:

```
example_notebooks/mdtf.natl_ocean.esnb.ipynb
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

The notebook needs three sets of files. On NCAR Casper they already exist at the paths below;
**off-NCAR users must copy them** (e.g. via Globus) and update the catalog + config paths.

**a) Reference / observation data** (point `OBS_DATA` at this directory):
- `obs_1x1.nc`
- `omip2.cycle1.1989_2018.0-200m.mld_sic_t_s_sigma.nc`
- `obs.maps_freq.sigma2.1982-2009_decomp_mean_1x1.nc`
- `obs_wmt_sigma2_1982-2009_spna_decomp_mean.nc`

On Casper: `/glade/work/taydral/MDTF/inputdata/obs_data/natl_ocean/`

**b) Model data** — for `data_source = "CESM_timeslice"`, the native POP files:
- `cesm_mdtfv3_timeslice.{TEMP,SALT,SHF,QFLUX,SFWF,UVEL,VVEL}.mon.nc`

On Casper: `/glade/campaign/cgd/amp/bundy/mdtf/cesm_mdtfv3_timeslice_public/ocn/mon/`
These paths are listed in the `path` column of
`diagnostics/natl_ocean/CESM_timeslice_ocean_monthly_001.csv`.

**c) fx (grid) data** (point `FX_DIR` at this directory) — CMIP6 CESM2 `areacello` and
`volcello` (gx1 grid, shared by both data sources):

On Casper: `/glade/collections/cmip/CMIP6/CMIP/NCAR/CESM2/historical/r1i1p1f1/Ofx/`

### Downloading via Globus (off-NCAR)

> **TODO (endpoint not yet created):** download the data sets above from the project Globus
> endpoint and place them in local directories of your choice.
> - Globus endpoint (model data): `TODO: <endpoint name / UUID>`
> - Globus endpoint (obs / reference data): `TODO: <endpoint name / UUID>`
> - Globus endpoint (fx data): `TODO: <endpoint name / UUID>`

### Updating catalog paths (off-NCAR)

The catalog files in `diagnostics/natl_ocean/` store **absolute** paths:
- the `path` column in `CESM_timeslice_ocean_monthly_001.csv` (and `CMIP_CESM_historical_001.csv`)
- the `catalog_file` field in the matching `.json`

After copying the model data locally, edit those `path` entries to point at your local copies
(or regenerate the catalog with MDTF's `catalog_builder` tool). The notebook resolves a
relative `catalog_file` against the catalog JSON's own directory, so a relative path works too.

---

## Step 4 — Configure paths

All site-specific paths live in one **USER CONFIGURATION** cell near the top of the notebook
(Section 1). Edit the `CHANGE ME` lines, or export the same variables in your shell before
launching Jupyter (the cell uses `setdefault`, so exported values win):

| Variable      | What to set it to |
|---------------|-------------------|
| `CODE_ROOT`   | your `MDTF-diagnostics` repo directory |
| `WORK_DIR`    | output / work directory (created if missing) |
| `OBS_DATA`    | the reference/obs directory from Step 3a |
| `FX_DIR`      | the fx directory from Step 3c |
| `PBS_ACCOUNT` | your PBS account (only used if you enable the dask cluster) |

---

## Step 5 — Launch and run

1. Start Jupyter:
   - **On Casper**: use JupyterHub (<https://jupyterhub.hpc.ucar.edu>) or a `qvscode`
     interactive session; a 1-year run is small enough to run single-process (no dask).
   - **Locally**: `conda activate esnb && jupyter lab`
2. Open `example_notebooks/mdtf.natl_ocean.esnb.ipynb` and select the **`Python (esnb)`**
   kernel (top-right kernel selector).
3. In Section 1, set `data_source` (`"CESM2"` or `"CESM_timeslice"`) and the date range.
   For a first test, use one year: `startdate, enddate = '1995-01-01', '1995-12-31'`
   (the timeslice data covers 1995–2015).
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

- MDTF-diagnostics docs: <https://mdtf-diagnostics.readthedocs.io/>
- Git setup guide: <https://mdtf-diagnostics.readthedocs.io/en/main/sphinx/dev_git_intro.html>
- POD detailed description: `diagnostics/natl_ocean/doc/natl_ocean.rst`
- ESNB: <https://github.com/jkrasting/esnb>
- intake-esm: <https://intake-esm.readthedocs.io/>

---

**POD authors / contributors:** Elizabeth Maroon, Stephen Yeager, Taydra Low, Brendan Myers, Feng Shu

**Last updated:** 2026-06-29
