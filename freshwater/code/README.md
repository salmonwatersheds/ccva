README
================
Stephanie Peacock
2025-12-24

## `freshwater/code`

Contains R code to summarise gridded hydrologic model output from the
Pacific Climate Impacts Consortium, available at , and calcualte the
exposure of 60 Conservation Units of salmon and steelhead in the Fraser
River basin.

### Files

**`1a_spatial-distribution-fraser.R`**: This code reads in a number of
spatial datasets related to salmon distribution in the **Fraser region**
and determines, for each CU, which PCIC grid cells should be used in
assessments of climate change exposure for freshwater life stages.

- Outputs:

  - `output/PCIC_incl.rds` - list with nrow = 60 CUs and ncol = 4
    freshwater life stages (adult_migration, spawning, eggs_alevin,
    fw_rearing), with each element is the unique PCIC grid cell ID that
    is considered habitat for that CU and sage.

  - `output/PCIC_incl_overall.csv` - simple vector of all unique PCIC
    grid cell ID that are habitat for any CU or stage, for subsetting
    hydrologic model output.

**`1b_process-PCIC_fraser.R`**: Code to load PCIC model output, trim to
spatial and temporal domain of interest output from
`1a_spatial-distribution-fraser.R`, calculate weekly rolling max or
average, and save output RDS file for use in
`2_calc-freshwater-exposure.R`

- Outputs:

  - `freshwater/data/processed-data/PCIC_CanESM2_rcp45_processed.rds`
    where `CanESM2` is one of six GCMs and `rcp45` can also be `rcp85`
    (12 total output files)

**`2_calc-freshwater-exposure.R`**: Code to summarize gridded hydrologic
model output from the PCIC and calculate projected changes in exposure
to stream temperature and flow outside of the optimal ranges for the
given species and life stage. Loops through each emissions scenario (n =
2), GCM (n = 6), CU (n = 60), and life stage (n = 4 in freshwater).

- Outputs:

  - `freshwater/output/freshwater_spat_fraser_rcp45_2024-12-19_shortSELmig.rds`:
    list with nrow = 60 CUs and ncol = 4 freshwater life stages
    (adult_migration, spawning, eggs_alevin, fw_rearing), with each
    element comtains another array of dim 6 (GCM), 2 (temp or flow), 4
    (period hist, early, mid, or late), **ncells** (number of grid cells
    relevant to the CU and stage).

  - `freshwater/output/freshwater_output_fraser_rcp45_2024-12-19_shortSELmig.rds`:
    an array with dim 6 (GCM), 60 (CUs), 4 (freshwater life stage), 2
    (temp or flow), 4 (period), 3 (median, lower, upper among all
    relevant grid cells).

  - Each of the above files are output for `rcp45` and `rcp85`.

  - Outputs are dated and denoted with `shortSELmig` if they take into
    account that lake-type sockeye can seek thermal refuge in lakes, and
    so the migration stage only includes the estimated time it takes to
    migrate from river entry to rearing lake (Table S4).

**`3_summarize-freshwater-exposure.R`**: Code to summarize projected
changes in exposure to stream temperature and flow outside of the
optimal ranges for the given species and life stage. Takes output from
`2_calc-freshwater-exposure.R` and tries to make sense of it! This code
is not very organized and final figures for paper are produced by code
in the `docs/` folder.
