README
================
Stephanie Peacock
2023-04-06

## Climate Change Vulnerability Assessments for Pacific Salmon and Steelhead

### Background

Climate Change Vulnerability Assessments (CCVAs) aim to quantify the
potential impact of climate change on certain valued attributes of a
particular system to inform the identification and prioritization of
conservation and mitigation actions. In this project, we are interested
in the potential impact of climate change on the viability of Pacific
salmon and steelhead Conservation Units (CUs), with the aim to
understand which CUs are likely to be most strongly affected and why.
For more information on the project, see the [Salmon Watersheds Program
website](https://salmonwatersheds.ca/project/ps13/).

### What’s in this repo

Our approach to CCVAs involves quantifying **exposure** to climate
changes as predicted by Global Climate Models, considering the unique
spatial distribution and life-history timing of each Conservation Unit.
We are considering climate changes throughout the salmon life cycle, and
the code and data here are organized into `freshwater` and `marine`
subfolders becuase the data sources and processing for these two
environments are different.

Within the `freshwater` and `marine` folders, there are R files in
`code` that read in climate output, organize and summarize it, and
calculate exposure for each CU, life stage, and climate attribute. **So
far, we have only completed these calculations for freshwater attributes
(stream temperature and low flow) for Fraser lake-type sockeye and are
working to expand to other species, regions, and the marine
environment.** This is a work in progress, so be kind if code isn’t
perfect and check back for updates!

Note that raw data files consist of Global Climate Model projections and
are large files (hundreds of GB). These are therefore not made available
on GitHub but can be downloaded from the relevant sources listed in the
table below or by contacting the project lead (Stephanie Peacock,
speacock at psf dot ca).

| Environment | Exposure attribute      | Definition                                                                                                                                                                                                  | Dataset                                                                                                                                                                                                                                                           |
|-------------|-------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Freshwater  | Stream temperature      | The proportion of days within each life stage that stream temperature is above the 90th percentile of historical temperatures experienced by that life stage                                                | [Pacific Climate Impacts Consortium gridded hydrologic model output](https://www.pacificclimate.org/data/gridded-hydrologic-model-output)                                                                                                                         |
| Freshwater  | Low flow                | The proportion of days within each life stage that flow is below 20% of the Mean Annual Discharge (MAD) for the grid cell                                                                                   | [Pacific Climate Impacts Consortium gridded hydrologic model output](https://www.pacificclimate.org/data/gridded-hydrologic-model-output)                                                                                                                         |
| Marine      | Sea surface temperature | The standardized change (i.e., z-score\*) from the historical period in sea surface temperature over the spatial domain of the early marine stage (CU-specific) or marine rearing stage (species-specific). | [NOAA’s Climate Change Web Portal](https://psl.noaa.gov/ipcc/ocn/ccwp.html)                                                                                                                                                                                       |
| Marine      | Sea surface salinity    | The standardized change (i.e., z-score\*) from the historical period in sea surface salinity over the spatial domain of the early marine stage (CU-specific) or marine rearing stage (species-specific).    | [NOAA’s Climate Change Web Portal](https://psl.noaa.gov/ipcc/ocn/ccwp.html)                                                                                                                                                                                       |
| Marine      | Sea surface pH          | The standardized change (i.e., z-score\*) from the historical period in sea surface pH over the spatial domain of the early marine stage (CU-specific) or marine rearing stage (species-specific).          | [NOAA’s Climate Change Web Portal](https://psl.noaa.gov/ipcc/ocn/ccwp.html)                                                                                                                                                                                       |
| Marine      | Sea level rise          | The predicted change in sea level along coastal migration routes of the early marine stage.                                                                                                                 | [James et al. (2021) Relative sea-level projections for Canada based on the IPCC Fifth Assessment Report and the NAD83v70VG national crustal velocity model](https://geoscan.nrcan.gc.ca/starweb/geoscan/servlet.starweb?path=geoscan/fulle.web&search1=R=327878) |

<!-- | Freshwater | Beneficial high flows | The change in return interval (in years) for an instantaneous high flow equal to a historical 10-year flood event. | [Pacific Climate Impacts Consortium gridded hydrologic model output](https://www.pacificclimate.org/data/gridded-hydrologic-model-output) | -->
<!-- | Freshwater | Damaging high flows | The change in the exceedance probability for an instantaneous high flow equal to a historical 100-year flood event. | [Pacific Climate Impacts Consortium gridded hydrologic model output](https://www.pacificclimate.org/data/gridded-hydrologic-model-output) | -->

### Draft outputs

The `docs` folder contains an R Shiny document that allows users to
explore draft results. Please note that this output is preliminary and
should not be over interpreted. The resulting Shiny document is online
at <https://salmonwatersheds.shinyapps.io/salmon-ccva/>.

## Acknowledgements

This work is funded by the BC Salmon Restoration and Innovation Fund. We
are extremely grateful to PSF’s Climate Science Advisory Committee, who
have provided input on our approach. We thank Markus Schnorbus at
[PCIC](www.pacificclimate.org) for providing preliminary hydrologic
model outputs so we can advance this work.
