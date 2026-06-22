# North Atlantic Ocean POD: ESNB Notebook User Guide

## Overview

The North Atlantic Ocean POD (`natl_ocean`) is a Process-Oriented Diagnostic for the MDTF framework that computes ocean heat and salt diagnostics for the subtropical-to-subpolar North Atlantic.

### What This POD Computes

- **Potential Density (σ₀)**: Stratification and density structure
- **Mixed Layer Depth (MLD)**: Surface mixing and water mass interaction
- **Water Mass Transformation (WMT)**: Transformation of water masses in density space (σ₂)
- **Atlantic Meridional Overturning Circulation (AMOC)**: Streamfunction in σ₂ coordinates
- **Bias Analysis**: Comparison against observations and model benchmarks

### Key Features

- Works with multiple data sources (CESM2, GFDL timeslice, other CMIP6 models)
- Two execution modes: **interactive** (manual) and **driver** (framework-driven)
- Two data loading approaches: **ESNB-based** and **raw intake-esm**
- Produces publication-ready diagnostic figures

---

## Part 1: Getting Started from Scratch

### Step 1: Clone the MDTF Repository

The natl_ocean POD is part of the MDTF-diagnostics package. Clone it from GitHub:

```bash
git clone https://github.com/NOAA-GFDL/MDTF-diagnostics.git
cd MDTF-diagnostics
```

### Step 2: Follow Git Setup Instructions

If this is your first time with MDTF development, read the git introduction:

https://mdtf-diagnostics.readthedocs.io/en/main/sphinx/dev_git_intro.html

Key steps:
- Configure your git user information
- Set up SSH keys (optional but recommended)
- Understand the branching model

### Step 3: Set Up Python Environment

The POD requires Python 3.12 with specific packages. Create a conda environment:

```bash
# Activate your miniconda installation
source /glade/u/home/taydral/miniconda3/etc/profile.d/conda.sh

# Create or activate the MDTF environment
conda activate _MDTF_python3_base

# Or create a fresh environment if needed:
# conda create -n mdtf-natl python=3.12 xarray dask intake-esm pandas matplotlib netCDF4 -c conda-forge
```

Required packages:
- xarray >= 0.20
- dask >= 2022.2.1
- intake-esm >= 2022.11.0
- pandas >= 1.4
- matplotlib >= 3.5
- netCDF4
- NumPy
- SciPy

### Step 4: Locate the Notebook

The natl_ocean ESNB notebook is at:

```
MDTF-diagnostics/example_notebooks/mdtf.natl_ocean_TL.ipynb
```

Two variants available:
- **mdtf.natl_ocean_TL.ipynb** — ESNB-based data loading (recommended for beginners)
- **mdtf.natl_ocean_intake_esm.ipynb** — Raw intake-esm API (for advanced users)

---

## Part 2: Running the Notebook

### Opening the Notebook

Start a Jupyter notebook server in the MDTF-diagnostics directory:

```bash
jupyter notebook example_notebooks/mdtf.natl_ocean_TL.ipynb
```

### Section 1: POD Settings and Configuration

This section defines the analysis parameters. Edit the cells to:

1. **Choose a data source** (cell: "Data Source Configuration")
   - Set `data_source = 'cesm2'` for CMIP6 historical data (default)
   - Set `data_source = 'gfdl_timeslice'` for GFDL timeslice model output

2. **Set the analysis period** (cell: "startdate, enddate")
   - Example: `startdate = '1985-01-01'`, `enddate = '2014-12-31'` (30 years)
   - Shorter periods (e.g., 1 year) work for testing but may not capture climate variability

3. **Configure computation options** (optional)
   - `regrid_resolution`: Target grid for regridding (default: 1.0°)
   - `zavg_depth`: Depth for averaging (default: 200m)

4. **Enable optional dask cluster** (optional, for HPC)
   - Uncomment the dask PBSCluster cell if running on Casper
   - Allows parallel computation of large datasets

### Section 2: Data Loading

This section loads data via intake-esm catalogs. **Run all cells in order without modification** for initial runs.

#### For CESM2 (default):

The notebook automatically uses the CMIP6 historical catalog. Data is loaded from:
```
/glade/collections/cmip/CMIP6/CMIP/NCAR/CESM2/historical/
```

#### For GFDL timeslice:

The notebook uses the GFDL timeslice catalog. Data is loaded from:
```
/glade/campaign/cgd/amp/bundy/mdtf/cesm_mdtfv3_timeslice_public/ocn/mon/
```

#### Important: Framework vs Interactive Mode

- **Interactive mode** (default): All data loading happens in the notebook cells
- **Driver mode** (used by MDTF framework): The framework handles data loading; skip Section 2

### Section 3: Diagnostics Computation (Steps 1-9)

This section is the core of the POD. It contains:

**Step 1: Preprocess coordinates**
- Standardize dimension names to POD conventions
- Compute layer thickness (dz)

**Step 2: Compute σ₀ and MLD**
- Potential density from temperature and salinity
- Mixed layer depth (depth where σ₀ exceeds surface by 0.03 kg/m³)

**Step 3: Depth-average fields**
- Average temperature, salinity, σ₀ in upper 200m
- Weighted by layer thickness

**Step 4: Regrid to regular lat-lon grid**
- Use xESMF for bilinear regridding to 1°×1° grid
- Enables comparison with observations on regular grid

**Step 5: Monthly climatology**
- Compute monthly means across the time period
- Tag with model name for later comparison

**Steps 6-7: Water Mass Transformation (WMT)**
- Prepare surface heat and salt fluxes
- Compute WMT rates in σ₂ (potential density referenced to 2000m) coordinates
- Diagnose transformation by region and density class

**Steps 8-9: AMOC calculation**
- Load velocity fields (U and V components)
- Calculate AMOC streamfunction in σ₂ coordinates
- Output: Atlantic Meridional Overturning Circulation magnitude

**To run**: Execute all cells in order. Most cells are independent and can be re-run if needed.

### Section 4: Results and Visualization

This section generates publication-ready figures:

1. **Bias maps**: Compare model vs. observations
   - Temperature, salinity, density, mixed layer depth
   - Shows regional biases and spatial patterns

2. **Scatter plots**: Temperature-salinity relationship
   - Reveals systematic biases in water mass characteristics

3. **WMT diagnostics**: Water mass transformation rates
   - By region (e.g., subpolar, subtropical)
   - By density class
   - By latitude (showing transformation hotspots)

4. **AMOC streamfunction**: Atlantic overturning circulation
   - Magnitude and latitude dependence
   - Comparison with observations (RAPID array at 26°N)

**Output**: All figures are saved as PNG files in the current directory.

---

## Part 3: Interpreting Results

### Understanding Potential Density (σ₀)

σ₀ = ρ(T, S, 0) - 1000, where:
- ρ = seawater density
- T = potential temperature
- S = salinity
- Reference pressure = 0 dbar (surface)

σ₀ > 27.8 indicates mode/intermediate waters

### Understanding WMT

Water Mass Transformation (WMT) rates (Sv/σ unit) show the volume of water crossing density surfaces due to surface forcing. High WMT indicates:
- Strong air-sea interaction
- Significant water mass modification
- Important regions for ocean circulation

### Understanding AMOC

AMOC streamfunction (Sv = 10⁶ m³/s):
- Positive values = clockwise circulation (northward heat transport)
- Peak typically occurs at 24-26°N (RAPID mooring latitude)
- Values: 15-20 Sv typical for realistic models

### Comparing Models

Biases indicate:
- Temperature bias → Incorrect heat content/stratification
- Salinity bias → Incorrect freshwater balance
- MLD bias → Mixing parameterization issues
- WMT bias → Surface forcing or mixing defects

---

## Part 4: Adapting the Notebook for Your POD

### Template Structure

This notebook serves as a template for creating new ocean PODs. The structure is:

1. **Section 1**: POD settings and knobs
   - User-configurable parameters
   - Data source selection
   - Computation options

2. **Section 2**: Data loading (framework-independent)
   - ESNB or intake-esm catalog access
   - Variable name mapping
   - Quality checks

3. **Section 3**: POD-specific diagnostics
   - Your physics calculations
   - Derived quantities
   - Preprocessing for output

4. **Section 4**: Visualization and output
   - Figure generation
   - Output saving
   - Result interpretation

### Steps to Create Your POD

1. **Copy this notebook** as a starting point:
   ```bash
   cp example_notebooks/mdtf.natl_ocean_TL.ipynb example_notebooks/mdtf.your_pod.ipynb
   ```

2. **Modify Section 1**: Update POD settings for your variables

3. **Modify Section 3**: Replace diagnostic computation with your physics

4. **Modify Section 4**: Change visualization to show your results

5. **Update settings.jsonc**: Declare variables your POD requires
   ```json
   {
     "varlist": {
       "your_var": {
         "frequency": "mon",
         "realm": "ocean",
         "dimensions": ["time", "lev", "nlat", "nlon"],
         "standard_name": "your_cmip_standard_name",
         "units": "units"
       }
     }
   }
   ```

6. **Test interactively**: Run the notebook to verify it works

7. **Test with framework**: Use MDTF framework to test driver mode

### Key Functions from POD_utils

The natl_ocean POD uses utility functions in `POD_utils.py`:

- `preprocess_coords(ds)` — Standardize coordinate names
- `compute_sigma0(da_t, da_s)` — Potential density
- `compute_mld(da_sig0)` — Mixed layer depth
- `compute_zavg(ds, var, dz, depth)` — Depth averaging
- `regrid(ds, target_grid)` — Regridding to regular grid
- `compute_wmt(ds, calc_type, density)` — Water mass transformation
- `calculate_moc(ds_t, ds_u, ds_v)` — AMOC streamfunction

For your POD, create similar utility functions for your own calculations.

### MDTF Framework Integration

When your POD is ready for framework integration:

1. Create `your_pod/your_driver.py` (or convert notebook to `.py`)
2. Create `your_pod/settings.jsonc` with variable declarations
3. Test with: `mdtf_framework.py your_pod/settings.jsonc`
4. Framework handles data preprocessing and provides variables as environment variables

Learn more: https://mdtf-diagnostics.readthedocs.io/

---

## Part 5: Troubleshooting

### Common Issues and Solutions

#### 1. ImportError: No module named 'esnb'

**Solution**: ESNB is only needed for the ESNB notebook variant. If using the raw intake-esm variant:
```python
import intake  # This should work
# Don't import: from esnb import ...
```

#### 2. Data not found in catalog

**Solution**: Verify the catalog path and variable names:
```python
cat = intake.open_esm_datastore(cat_file)
print(cat.df)  # Shows all available variables
```

#### 3. Dask cluster fails

**Solution**: Fall back to single-process mode by commenting out dask setup:
```python
# client = dask_jobqueue.PBSCluster(...)
# Computations will use single CPU instead
```

#### 4. Memory errors with large datasets

**Solutions**:
- Reduce time period (e.g., 1 year instead of 30 years)
- Use higher regridding resolution (coarser grid)
- Enable dask for parallel computation
- Process one variable at a time

#### 5. Coordinate mismatch errors

**Solution**: The notebook includes `preprocess_coords()` to handle different data sources:
```python
ds_target = POD_utils.preprocess_coords(ds_target)
# This standardizes dimension names (lat→nlat, lon→nlon, etc.)
```

---

## Part 6: References and Further Reading

### Official Documentation
- MDTF-diagnostics: https://mdtf-diagnostics.readthedocs.io/
- Git setup guide: https://mdtf-diagnostics.readthedocs.io/en/main/sphinx/dev_git_intro.html
- POD development guide: https://mdtf-diagnostics.readthedocs.io/en/main/sphinx/dev_guidelines.html

### Data and Tools
- intake-esm: https://intake-esm.readthedocs.io/
- ESNB (Enhanced Jupyter Notebooks): https://github.com/NOAA-GFDL/esnb
- xarray: https://xarray.pydata.org/
- dask: https://dask.org/

### Scientific References
- CMIP6 standard names: https://cmip.llnl.gov/
- RAPID AMOC array: https://www.rapid.ac.uk/
- Water mass transformation: Speer and Tziperman (1992), J. Phys. Oceanogr.
- AMOC in density coordinates: Gent (2016), J. Phys. Oceanogr.

---

## Appendix: Data Sources

### CESM2 CMIP6 Historical (Default)

**Location**: `/glade/collections/cmip/CMIP6/CMIP/NCAR/CESM2/historical/`

**Catalog**: CMIP_CESM_historical_001.json

**Variables**: thetao, so, uo, vo, tos, hfds, vsf (+ static: areacello, volcello)

**Time range**: 1850-2014

**Resolution**: Native CESM2 (~0.9° × 1.25°)

### GFDL Timeslice Data (Alternative)

**Location**: `/glade/campaign/cgd/amp/bundy/mdtf/cesm_mdtfv3_timeslice_public/ocn/mon/`

**Catalog**: GFDL_timeslice_ocean_monthly_001.json

**Variables** (catalog `variable_id` → native POP name): thetao→TEMP, so→SALT, uo→UVEL,
vo→VVEL, hfds→QFLUX, vsf→SFWF. The files store the *native POP* variable names; the catalog
indexes them under CMIP `variable_id` so the same queries work, but the loaded datasets carry
POP variable/coordinate names.

**Format**: Aggregated monthly files (1 file per variable), 1995–2015 (240 months)

**Resolution**: Native CESM/POP gx1 displaced-pole grid (384×320, `grid_label=gn`); depth
`z_t` in **centimeters**; horizontal coords `TLONG/TLAT` (T-grid) and `ULONG/ULAT` (U-grid).
This is *not* a regridded 1° product.

> **Status note**: The GFDL catalog opens and resolves all six variables, and selecting
> `data_source = 'gfdl_timeslice'` loads them. However, because the data are raw POP output,
> the GFDL branch still needs a harmonization step before the full diagnostic pipeline
> (σ₀/MLD/zavg/WMT/AMOC) will run: rename data vars to CMIP names, rename `z_t→lev`,
> `TLONG/TLAT→lon/lat`, and construct `lev_bnds` from `z_t`. End-to-end GFDL validation is
> pending. The CESM2 path is fully validated.

---

**Last Updated**: June 2026

**Notebook Authors**: Taydra Lowenstein, Liz Hertz (GFDL/NCAR)

**Questions?** Contact the MDTF development team or create an issue on GitHub:
https://github.com/NOAA-GFDL/MDTF-diagnostics/issues
