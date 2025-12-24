README
================
Stephanie Peacock
2025-12-24

## `marine/code`

Contains R code to summarise raw GCM output from CMIP5 climate change
model projections, provided by NOAA at
<https://psl.noaa.gov/ipcc/ocn/ccwp.html>, and calculate the exposure of
60 Conservation Units of salmon and steelhead in the Fraser River basin.

### Files

**`1a_subset-NOAA-grid.R`**: This code reads in a number of spatial
datasets related to determine which grid points for the coarse GCM model
output provided by NOAA should be included for the Fraser early marine
and marine rearing stages.

- Outputs:

  - `output/incl_NOAA_marRear.csv` - simple vector of all unique NOAA
    grid cell IDs that are habitat for salmon in the **marine rearing**
    stage (assumed the same for all CUs).
  - `output/incl_NOAA_RM.csv` - simple vector of all unique NOAA grid
    cell IDs that are habitat for salmon in the **early marine** stage
    (assumed the same for all CUs).

**`2_calc-marine-exposure.R`**: Code to summarize gridded climate change
projections for sea surface temp and salinity from CMIP5 GCMs and
calculate projected changes in exposure to temperature and salinity
outside of the optimal range for each species(for temperature) and
historical norms in each grid cell (for salinity).

- Outputs:

  - `marine/output/marine_spat_fraser_rcp45_2024-02-06.rds`: list with
    nrow = 60 CUs and ncol = 2 marine life stages (early_marine,
    marine_rearing), with each element contains another array of dim 6
    (GCM), 2 (temp or salinity), 4 (period hist, early, mid, or late),
    **ncells** (number of grid cells relevant to the stage (same for all
    CUs in marine stages)).

  - `marine/output/marine_output_fraser_rcp45_2024-02-06.rds`: array
    with dim 6 (GCMs), 60 (CUs), 2 (early_marine or marine_rearing), 2
    (temp or salinity), 4 (period hist, early, mid, or late), 3 (median,
    lower, upper among all grid cells).
