README
================
Stephanie Peacock
2023-04-06

## Climate Change Vulnerability Assessments for Pacific Salmon and Steelhead

Climate Change Vulnerability Assessments (CCVAs) aim to quantify the
potential impact of climate change on certain valued attributes of a
particular system to inform the identification and prioritization of
conservation and mitigation actions. In this project, we are interested
in the potential impact of climate change on the viability of Pacific
salmon and steelhead Conservation Units (CUs), with the aim to
understand which CUs are likely to be most strongly affected and why.

Our approach to CCVAs involves quantifying **exposure** to climate
changes as predicted by Global Climate Models, considering the unique
spatial distribution and life-history timing of each Conservation Unit.
The exposure attributes that we consider are listed in the following
table. Each attribute is summarized by CU and life stage by the
associated .R file in `code`, with resulting exposure saved in
`outputs`. Raw data files consist of Global Climate Model projections
and are large files (hundreds of GB). These are therefore not made
available on GitHub but can be downloaded from the relevant sources
listed in the table below or by contacting the project lead (Stephanie
Peacock, speacock at psf dot ca).

| Environment | Exposure attribute            | Definition                                                                                                                                                                                                  | Dataset                                                                                                                                                                                                                                                           |
|-------------|-------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Freshwater  | Stream temperature - optimal  | The change from the historical period in the number of days within the relevant life-history period when average stream temperature is outside the optimal temperature range for that life stage            | [Pacific Climate Impacts Consortium gridded hydrologic model output](https://www.pacificclimate.org/data/gridded-hydrologic-model-output)                                                                                                                         |
| Freshwater  | Stream temperature - critical | The change from the historical period in the number of days within the relevant life-history period when average stream temperature is outside the optimal temperature range for that life stage            | [Pacific Climate Impacts Consortium gridded hydrologic model output](https://www.pacificclimate.org/data/gridded-hydrologic-model-output)                                                                                                                         |
| Freshwater  | Low flow                      | The number of days within the relevant life-history period when the daily flow is below 20% historical mean annual discharge (MAD).                                                                         | [Pacific Climate Impacts Consortium gridded hydrologic model output](https://www.pacificclimate.org/data/gridded-hydrologic-model-output)                                                                                                                         |
| Marine      | Sea surface temperature       | The standardized change (i.e., z-score\*) from the historical period in sea surface temperature over the spatial domain of the early marine stage (CU-specific) or marine rearing stage (species-specific). | [NOAA’s Climate Change Web Portal](https://psl.noaa.gov/ipcc/ocn/ccwp.html)                                                                                                                                                                                       |
| Marine      | Sea surface salinity          | The standardized change (i.e., z-score\*) from the historical period in sea surface salinity over the spatial domain of the early marine stage (CU-specific) or marine rearing stage (species-specific).    | [NOAA’s Climate Change Web Portal](https://psl.noaa.gov/ipcc/ocn/ccwp.html)                                                                                                                                                                                       |
| Marine      | Sea surface pH                | The standardized change (i.e., z-score\*) from the historical period in sea surface pH over the spatial domain of the early marine stage (CU-specific) or marine rearing stage (species-specific).          | [NOAA’s Climate Change Web Portal](https://psl.noaa.gov/ipcc/ocn/ccwp.html)                                                                                                                                                                                       |
| Marine      | Sea level rise                | The predicted change in sea level along coastal migration routes of the early marine stage.                                                                                                                 | [James et al. (2021) Relative sea-level projections for Canada based on the IPCC Fifth Assessment Report and the NAD83v70VG national crustal velocity model](https://geoscan.nrcan.gc.ca/starweb/geoscan/servlet.starweb?path=geoscan/fulle.web&search1=R=327878) |

<!-- | Freshwater | Beneficial high flows | The change in return interval (in years) for an instantaneous high flow equal to a historical 10-year flood event. | [Pacific Climate Impacts Consortium gridded hydrologic model output](https://www.pacificclimate.org/data/gridded-hydrologic-model-output) | -->
<!-- | Freshwater | Damaging high flows | The change in the exceedance probability for an instantaneous high flow equal to a historical 100-year flood event. | [Pacific Climate Impacts Consortium gridded hydrologic model output](https://www.pacificclimate.org/data/gridded-hydrologic-model-output) | -->
