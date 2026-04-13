# MDTF-diagnostics

## Project Overview
The MDTF-diagnostics package is a portable framework for running process-oriented diagnostics (PODs) on weather and climate model data. Each POD targets a specific physical process or emergent behavior to determine how well models represent the process, ensure models produce the right answers for the right reasons, and identify gaps in understanding. PODs generate diagnostic figures viewable as HTML in a web browser.

## Resources
- github.com/NOAA-GFDL/MDTF-diagnostics
- mdtf-diagnostics.readthedocs.io/en/main/

## Environment
- **Machine**: NCAR Casper HPC (PBS job scheduler)
- **Python**: 3.12 via miniconda3 (`/glade/u/home/taydral/miniconda3/`)
- **Primary conda env**: `_MDTF_python3_base`
- **Job scheduler**: PBS (`qsub`, `qstat`, `qdel`)
- **Key libraries**: xarray, dask, pandas, numpy, matplotlib, yaml, intake-esm
- **ESNB framework**: Notebook-based diagnostic execution via `esnb` package

## Project Structure
- `example_notebooks/` — Jupyter notebooks for developing and testing PODs (e.g., blocking diagnostic)
- `natl_pod/` — North Atlantic POD development

## Data Paths
- Model data: `/glade/campaign/cesm/` and `/glade/campaign/cgd/`
- Obs data: `/glade/work/bundy/mdtf/inputdata/obs_data/` and `/glade/work/rneale/data/`
- Output: `/glade/work/bundy/mdtf/outdir/`
- MDTF catalogs: `/glade/u/home/bundy/diag/mdtf/catalogs/`

## Conventions
- POD settings are defined in `settings.jsonc` files or inline dicts in notebooks
- Variable naming follows CMIP/CESM conventions (e.g., `zg` for geopotential height)
- Use `blocking_utils` and `blocking_figs` modules for blocking diagnostic analysis and plotting
- Notebooks support two modes: `interactive` (Jupyter) and `driver` (command line)
- PBS jobs use account codes like `UWIS0041`

## HPC Notes
- Use `qsub` to submit batch jobs, `qstat -u $USER` to check status
- Dask clusters via `dask_jobqueue.PBSCluster` for parallel computation
- Large datasets live on `/glade/campaign/` (read-only) and `/glade/work/` (read-write)