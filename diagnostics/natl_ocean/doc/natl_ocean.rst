North Atlantic Ocean Diagnostic Documentation
=============================================

Last update: 6/29/2026

This POD evaluates how well a model represents the subtropical-to-subpolar North
Atlantic ocean state, comparing the model against a 1°×1° observation-based reference.
It diagnoses upper-ocean temperature and salinity, mixed layer depth, potential density,
the surface heat and freshwater forcing, and the density-space circulation that the
forcing drives. Specifically it produces:

- Sea-surface and upper-200 m temperature/salinity climatologies and their biases,
- Mixed layer depth (MLD) and potential density (σ₀) climatologies,
- Surface heat and freshwater flux fields,
- Surface-forced water-mass transformation (WMT) in σ₂ coordinates, by region,
- The Atlantic meridional overturning circulation (AMOC) streamfunction in σ₂ coordinates,
- A WMT + AMOC-at-45°N compensation diagnostic relating the surface-forced transformation
  to the overturning.

The diagnostic is data-source agnostic: it runs on CMIP-conformant (CMORized) output
(e.g. CESM2 historical, loaded through the ESNB framework) and on native model output
that is harmonized to CMIP variable/coordinate conventions in the example notebook
(e.g. CESM "timeslice" POP output, for which ``tos`` is derived from surface ``thetao``).
For native CESM/POP2 output only, the notebook applies the standard frazil surface-flux
recipe — ``hfds = SHF + QFLUX`` and ``wfo = SFWF - QFLUX/(latent_heat_fusion/1e4)`` — so
that the surface heat and freshwater forcing match CMIP conventions; CMORized inputs are
already corrected and skip this step.

Version & Contact info
----------------------

- Version/revision information: version 1 (6/29/2026)
- PIs: Liz Maroon (University of Wisconsin, emaroon@wisc.edu) and Steve Yeager (NSF National Center for Atmospheric Research, yeager@ucar.edu)
- Developer/point of contact: Liz Maroon (University of Wisconsin, emaroon@wisc.edu)
- Other contributors: Taydra Low, Brendan Myers, Teagan King

Open source copyright agreement
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The MDTF framework is distributed under the LGPLv3 license (see LICENSE.txt).
Unless you've distributed your script elsewhere, you don't need to change this.

Functionality
-------------

1. Preprocessing & Derived Variable Computation:
     Prepare the model output and compute derived physical variables using POD_utils.py.
     Key Functions:
       preprocess_coords(ds): Harmonizes coordinate/dimension names (nlat/nlon → y/x, latitude/longitude → lat/lon, etc.).
       check_depth_units(ds): Ensures depth (lev) is in meters, converting from cm if needed.
       compute_sigma0(thetao, so): Calculates potential density anomaly (σ₀).
       compute_mld(sigma0): Computes Mixed Layer Depth based on a density threshold.
       compute_zavg(ds, var): Calculates thickness-weighted mean of a variable over the upper 200 m.

2. Regridding:
     Interpolate all model and observational fields to a consistent 1×1 lat-lon grid using POD_utils.py.
     Key Functions:
       regrid(ds, method='bilinear'): Uses xESMF for bilinear regridding, with NaN handling and global domain definition.

3. Time Averaging (Climatology)
     Convert the time series into a monthly climatology over the analysis period.
     Note that it is recommended to use climo years that match the OMIP benchmark file, but this is not required.

4. Bias Calculation Against Observations
     Compute model bias relative to observations and calculate error statistics using POD_utils.py.
     Key Functions:
       plot_preproc(...): Extracts spatial and temporal slices, sets up metadata.
       error_stats(...): Computes area-weighted RMSE (map region) and mean bias (focus region).

5. Surface-Forced Water-Mass Transformation (WMT)
     Compute the transformation of water across σ₂ density classes driven by surface heat
     and freshwater fluxes, integrated over the subpolar North Atlantic, using POD_utils.py.
     Key Functions:
       compute_wmt(...): Surface-forced transformation rate as a function of σ₂, decomposed
         into heat and freshwater contributions; produces both regional WMT lines and σ₂
         density-flux maps.

6. AMOC Streamfunction in σ₂ Coordinates
     Compute the Atlantic overturning streamfunction in density space from the velocity
     (or mass-transport) fields, using POD_utils.py (AMOC code adapted from Sub2Sub's
     moc_funcs.py).
     Key Functions:
       calculate_moc(...): Integrates the meridional transport in σ₂ space. Accepts either
         velocity inputs (uo/vo, use_currents=True) or mass-transport inputs (umo/vmo,
         use_currents=False), handling the U-grid vs T-grid distinction.

7. Plotting
     Generate spatial and statistical visualizations of model performance using POD_utils.py.
     Key Functions:
       SpatialPlot_climo_bias(...): Two-panel maps for each variable (bias + bias rank).
       ScatterPlot_Error(...): Scatter plot comparing regional bias statistics between two variables.
       wmt_plot_byregion(...) / wmt_plot_maps(...): WMT lines vs OMIP benchmarks and σ₂ density-flux maps.

8. Optional: Multi-Cycle OMIP Time Handling
     Detect and restructure repeating OMIP forcing cycles.
     Key Functions:
       forcing_cycles(expid, nt): Identifies the number and span of OMIP forcing cycles.

Required programming language and libraries
-------------------------------------------

The North Atlantic Ocean Diagnostic recommends python (3.10 or later) because we
use xarray. Xarray, matplotlib, os, yaml, intake, numpy, xesmf, xskillscore,
scipy, gsw_xarray, numba, cftime, pandas, and cartopy are also required.

Required model output variables
-------------------------------
tos        Sea surface temperature       degrees Celsius   3D: time × lat × lon  (may be derived from surface thetao)
thetao     Potential temperature         degrees Celsius   4D: time × depth × lat × lon
so         Salinity                      PSU               4D: time × depth × lat × lon
hfds       Surface downward heat flux    W m-2             3D: time × lat × lon
wfo / vsf  Surface freshwater / virtual  kg m-2 s-1        3D: time × lat × lon  (at least one required)
           salt flux
uo         Sea water x velocity          m s-1             4D: time × depth × lat × lon
vo         Sea water y velocity          m s-1             4D: time × depth × lat × lon
umo, vmo   Ocean mass x/y transport      kg s-1            4D  (optional; enables the mass-transport AMOC path)
areacello  Ocean grid cell area          m2                2D: lat × lon  (fx)
volcello   Ocean grid cell volume        m3                3D: depth × lat × lon  (fx)
lev        Depth level                   m or cm           1D
lon, lat   Longitude and latitude        degrees           1D or 2D grid coordinates
time       Time dimension                datetime          1D

References
----------

.. _ref-Maloney:

1. E. D. Maloney et al. (2019): Process-Oriented Evaluation of Climate and
   Weather Forecasting Models. *BAMS*, **100** (9), 1665–1686,
   `doi:10.1175/BAMS-D-18-0042.1 <https://doi.org/10.1175/BAMS-D-18-0042.1>`__.

.. TODO: The science papers describing this diagnostic are in preparation
   (Maroon, Yeager, Low, and collaborators). Add the full citation(s) and DOI(s)
   here once they are published.

More about this diagnostic
--------------------------

Surface-forced water-mass transformation (WMT) quantifies the rate at which surface
buoyancy fluxes — heat and freshwater — convert water from one density class to another.
In the subpolar North Atlantic, surface cooling transforms light water into the dense
water that feeds the lower limb of the AMOC. Comparing the surface-forced transformation
in σ₂ space against the AMOC streamfunction in the same density coordinate (the WMT +
AMOC-at-45°N compensation diagnostic) tests whether a model's overturning is consistent
with the water-mass transformation its surface forcing implies. Persistent biases in
upper-ocean temperature/salinity, MLD, or the surface fluxes therefore propagate into
biases in both the transformation and the overturning, which this POD makes visible by
plotting them side by side against the OMIP benchmark and the observation-based reference.
