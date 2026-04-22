North Atlantic Ocean Diagnostic Documentation
=============================================

Last update: 6/3/2025

TODO: synopsis of your diagnostic here (like an abstract). 

TODO: test by copying the source file into the online editor 
at `https://livesphinx.herokuapp.com/ <https://livesphinx.herokuapp.com/>`__ and 
experiment.

Version & Contact info
----------------------

- Version/revision information: version 1 (6/3/2025)
- PIs: Liz Maroon (University of Wisconsin, emaroon@wisc.edu) and Steve Yeager (NSF National Center for Atmospheric Resarch, yeager@ucar.edu)
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
       compute_sigma0(thetao, so): Calculates potential density anomaly (σ₀).
       compute_mld(sigma0): Computes Mixed Layer Depth based on a density threshold.
       compute_zavg(ds, var): Calculates thickness-weighted mean of a variable over the upper 200 m.
       check_depth_units(ds): Ensures depth units are in meters.

2. Regridding:
     Interpolate all model and observational fields to a consistent 1×1 lat-lon grid using POD_utils.py.
     Key Functions:
       regrid(ds, method='bilinear'): Uses xESMF for bilinear regridding, with NaN handling and global domain definition.

3. Time Averaging (Climatology)
     Convert the time series into a monthly climatology over 1989–2018 using user-controlled driver script.
     Note that it is recommended to use climo_years from the model experiment that match the OMIP file, but not necessary.

4. Bias Calculation Against Observations
     Compute model bias relative to observations and calculate error statistics using POD_utils.py.
     Key Functions:
       plot_preproc(...): Extracts spatial and temporal slices, sets up metadata.
       error_stats(...): Computes area-weighted RMSE (map region) and mean bias (focus region).

5. Plotting
     Generate spatial and statistical visualizations of model performance using POD_utils.py.
     Key Functions:
       SpatialPlot_climo_bias(...): Produces two-panel maps for each variable (bias + bias rank).
         SpatialBias_panel(...): Shows spatial bias and regional stats.
         SpatialRank_panel(...): Shows where the target model ranks in bias among OMIP models.
       ScatterPlot_Error(...): Scatter plot comparing regional bias statistics between two variables.
         Uses Scatter_panel(...) to visualize model-model comparisons.

6. Optional: Multi-Cycle OMIP Time Handling
     Detect and restructure repeating OMIP forcing cycles.
     Key Functions:
       forcing_cycles(expid, nt): Identifies number and span of OMIP cycles.
       reorg_by_cycle(ds, nt, ncyc, yearrange): Restructures the dataset to expose cycles on a new dimension.

Required programming language and libraries
-------------------------------------------

The North Atlantic Ocean Diagnostic recommends python (3.10 or later) because we
use xarray. Xarray, matplotlib, os, yaml, intake, numpy, xesmf, xskillscore,
scipy, gsw_xarray, numba, cftime, and cartopy are also required.

Required model output variables
-------------------------------
thetao    Potential temperature   degrees Celsius  3D: time × depth × lat × lon
so        Salinity                PSU              3D: time × depth × lat × lon
lev       Depth level             m or cm          2D
lev_bnds  Depth level bounds      m or cm          2D: lev × bnds (bnds = 2)
lon, lat  Longitude and latitude  degrees          1D or 2D grid coordinates
time      Time dimension          datetime         1D

References
----------

TODO: Here you should cite the journal articles providing the scientific basis for 
your diagnostic. To keep the documentation format used in version 2.0 of
the framework, we list references "manually" with the following command:

.. Note this syntax, which sets the "anchor" for the hyperlink: two periods, one
   space, one underscore, the reference tag, and a colon, then a blank line.

.. code-block:: restructuredtext

   .. _ref-Maloney: 

   1. E. D. Maloney et al. (2019): Process-Oriented Evaluation of Climate 
   and Weather Forecasting Models. *BAMS*, **100** (9), 1665–1686, 
   `doi:10.1175/BAMS-D-18-0042.1 <https://doi.org/10.1175/BAMS-D-18-0042.1>`__.

which produces

.. _ref-Maloney: 
   
1. E. D. Maloney et al. (2019): Process-Oriented Evaluation of Climate and 
Weather Forecasting Models. *BAMS*, **100** (9), 1665–1686, 
`doi:10.1175/BAMS-D-18-0042.1 <https://doi.org/10.1175/BAMS-D-18-0042.1>`__.

which can be cited in text as ``:ref:`a hyperlink <reference tag>```, which 
gives :ref:`a hyperlink <ref-Maloney>` to the location of the reference on the 
page. Because references are split between this section and the following "More 
about this diagnostic" section, unfortunately you'll have to number references 
manually.

We don't enforce any particular bibliographic style, but please provide a 
hyperlink to the article's DOI for ease of online access. Hyperlinks are written
as ```link text <URL>`__`` (text and url enclosed in backticks, followed by two 
underscores).

More about this diagnostic
--------------------------

In this section, you can go into more detail on the science behind your 
diagnostic, for example, by copying in relevant text articles you've written. 
It's especially helpful if you're able to teach users how to use 
your diagnostic's output, by showing how to interpret example plots.

Instead of doing that here, we provide more examples of RestructuredText
syntax that you can customize as needed.

As mentioned above, we recommend the online editor at `https://livesphinx.herokuapp.com/ 
<https://livesphinx.herokuapp.com/>`__, which gives immediate feedback and has
support for sphinx-specific commands.


Links to external sites
^^^^^^^^^^^^^^^^^^^^^^^

URLs written out in the text are linked automatically: https://ncar.ucar.edu/. 

To use custom text for the link, use the syntax 
```link text <https://www.noaa.gov/>`__`` (text and url enclosed in backticks, 
followed by two underscores). This produces `link text <https://www.noaa.gov/>`__.

More references and citations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here's another reference:

.. code-block:: restructuredtext

   .. _ref-Charney: 

   2. Charney, Jule; Fjørtoft, Ragnar; von Neumann, John (1950). Numerical 
   Integration of the Barotropic Vorticity Equation. *Tellus* **2** (4) 237–254, 
   `doi:10.3402/tellusa.v2i4.8607 <https://doi.org/10.3402/tellusa.v2i4.8607>`__.

.. _ref-Charney: 

2. Charney, Jule; Fjørtoft, Ragnar; von Neumann, John (1950). Numerical 
Integration of the Barotropic Vorticity Equation. *Tellus* **2** (4) 237–254, 
`doi:10.3402/tellusa.v2i4.8607 <https://doi.org/10.3402/tellusa.v2i4.8607>`__.

Here's an example of citing these references:

.. code-block:: restructuredtext

   :ref:`Maloney et. al., 2019 <ref-Maloney>`, 
   :ref:`Charney, Fjørtoft and von Neumann, 1950 <ref-Charney>`

produces :ref:`Maloney et. al., 2019 <ref-Maloney>`, 
:ref:`Charney, Fjørtoft and von Neumann, 1950 <ref-Charney>`.

Figures
^^^^^^^

Images **must** be provided in either .png or .jpeg formats in order to be 
displayed properly in both the html and pdf output.

Here's the syntax for including a figure in the document:

.. code-block:: restructuredtext

   .. _my-figure-tag: [only needed for linking to figures]

   .. figure:: [path to image file, relative to the source.rst file]
      :align: left
      :width: 75 % [these both need to be indented by three spaces]

      Paragraphs or other text following the figure that are indented by three
      spaces are treated as a caption/legend, eg:

      - red line: a Gaussian
      - blue line: another Gaussian

which produces

.. _my-figure-tag:

.. figure:: gaussians.jpg
   :align: left
   :width: 75 %

   Paragraphs or other text following the figure that are indented by three
   spaces are treated as a caption/legend, eg:

   - blue line: a Gaussian
   - orange line: another Gaussian

The tag lets you refer to figures in the text, e.g. 
``:ref:`Figure 1 <my-figure-tag>``` → :ref:`Figure 1 <my-figure-tag>`.
