###############################################################################
# Code to extract which grid cells of the PCIC gridded hydrologic model output 
# from the Pacific Climate Impacts Consortium are used when summarizing 
# freshwater exposure for each CU and life stage.
# 
# Contact: Stephanie Peacock <speacock@psf.ca>
# Date: Sept 15, 2023
###############################################################################

# ** Currently works for Fraser SEL **

#------------------------------------------------------------------------------
# List of CUs in the PSE
#------------------------------------------------------------------------------

cu_list <- read.csv("data/cu_list.csv") %>% 
  subset(Species == "Lake sockeye" & Region == "Fraser" & Notes != "Extinct")

n.CUs <- length(unique(cu_list$Conservation.Unit))

#------------------------------------------------------------------------------
# Load full PCIC grid
#------------------------------------------------------------------------------

grid_points <- read.csv("output/PCIC-grid-points_Fraser.csv") %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4269)

#------------------------------------------------------------------------------
# Spawning Zone of Influence: Fraser SEL
#------------------------------------------------------------------------------

# Downloaded from Data Libary Sept 18, 2023
# https://data.salmonwatersheds.ca/result?datasetid=112
spawn_zoi <- st_read(dsn = "data/spatial/ZOI/se_spwng_zoi/se_spwng_zoi.shp") %>% st_transform(crs = 4269)

#------------------------------------------------------------------------------
# Adult migration lines: Fraser SEL
#------------------------------------------------------------------------------
mig_paths <- st_read(dsn = "data/spatial/fw-migration/spawn_timing_migration_paths_wCU_fraser.shp") %>% st_transform(crs = 4269)

#------------------------------------------------------------------------------
# Conservation Unit Boundaries: SEL
#------------------------------------------------------------------------------

# DFO shapefiels from Open Data Canada
# shp_path <- "data/spatial/CU-boundaries/Lake_Type_Sockeye_Salmon_CU_Boundary/SEL_CU_BOUNDARY_En"
# cu_boundary <- st_read(dsn = paste0(shp_path, ".shp"), promote_to_multi = FALSE)

# PSF version of shapefiles
shp_path <- "data/spatial/CU-boundaries/PSF_CUs_updatedDec2021"
cu_boundary <- st_read(dsn = paste0(shp_path, ".shp"), promote_to_multi = FALSE) %>%
  subset(regionname == "Fraser" & species %in% c("Lake sockeye", "Sockeye-Lake"))

#-----------------------------------------------------------------------------
# Define freshwater stages
#-----------------------------------------------------------------------------

stages <- c("adult_migration", "spawning", "eggs_alevin", "fw_rearing")

#-----------------------------------------------------------------------------
# Read in example of PCIC data to exclude cells that are NA
#-----------------------------------------------------------------------------
z <- nc_open(paste0("data/raw-data/PCIC/hydro_model_out/fraser/", "waterTemp", "_day_dynWat-VICGL_", "CanESM2", "_rcp45_", "r1i1p1", "_19450101-20991231_fraser.nc"))
varExtracted0 <-  ncvar_get(z, "waterTemperature")

# Spatial variables
lon <- ncvar_get(z, "lon")
lat <- ncvar_get(z, "lat")

# Create time variable: all output runs from 19450101 to 20991231
date <- as.Date(c(0:56612), origin = "1945-01-01")

# Reduce dimensionality so that grid point[i, ] corresponds to Ts[i, ]
varExtracted1 <- varExtracted0
dim(varExtracted1) <- c(length(lon)*length(lat), length(date))

# Historical period is 1970-2000, so no need to keep data pre-1970
date_keep <- which(date >= as.Date("1970-01-01") - 7)
Ts <- varExtracted1[, date_keep] 

rm(varExtracted1, varExtracted0)

#------------------------------------------------------------------------------
# Dates since 1970-01-01 minus 7 days (for averaging)
#------------------------------------------------------------------------------

# Historical period is 1970-2000, so no need to keep data pre-1970
date_keep <- which(date >= as.Date("1970-01-01") - 7)
#-----------------------------------------------------------------------------
# Store grid polys for each life stage and CU
#-----------------------------------------------------------------------------

# List of total c
incl.stages <- list(); length(incl.stages) <- n.CUs*4; dim(incl.stages) <- c(n.CUs, 4)

for(i in 1:n.CUs){
  
  cuid <- cu_list$CUID[i]
  
  # Subset spatial data for selected CU
  
  if(length(which(cu_boundary$FULL_CU_IN == cu_list$Full.CU.Index[i])) == 0){
    stop(paste0("No CU boundary for ", cu_list$Conservation.Unit[i]))
  } else {
    cu_boundary.i <- cu_boundary[which(cu_boundary$FULL_CU_IN == cu_list$Full.CU.Index[i]),]
  }
  
  zoi.i <- spawn_zoi[which(spawn_zoi$cuid == cuid), ]
  
  mig_paths.i <- mig_paths[which(mig_paths$cuid == cuid), ]
  
  for(j in 1:4){ # For each freshwater life stage
    
    # For adult migration and freshwater rearing, include downstream migration route
    if(stages[j] %in% c("adult_migration", "fw_rearing")){
      intrscts <- st_is_within_distance(grid_points, mig_paths.i, 3000, sparse = FALSE)  
      incl.stages[[i,j]] <- which(c(intrscts[, 1]) == TRUE)
      # There may be more than one line segment
      if(dim(intrscts)[2] > 1){
        for(k in 2:dim(intrscts)[2]){
          incl.stages[[i,j]] <- unique(c(incl.stages[[i,j]], which(c(intrscts[, k]) == TRUE)))
        }
      }
    } 
    
    # For spawning, eggs_alevin, use spawning ZOI
    if(stages[j] %in% c("spawning", "eggs_alevin")){
      intrscts <- st_intersects(grid_points, zoi.i, sparse = FALSE)
      incl.stages[[i,j]] <- which(c(intrscts) == TRUE)
    }
    
    # For freshwater rearing, include rearing lake (CU boundary)
    if(stages[j] == "fw_rearing"){
      
      # CU boundary is a multipolygon and was throwing error: 
      p <- st_make_valid(cu_boundary.i)
      
      intrscts <- st_intersects(grid_points, p, sparse = FALSE)
      incl.stages[[i,j]] <- unique(c(incl.stages[[i,j]], which(c(intrscts) == TRUE)))
    } 
    
    # Exclude cells that have no data for any years
    incl.stages[[i,j]] <- incl.stages[[i,j]][apply(Ts[incl.stages[[i,j]], ], 1, function(x) sum(!is.na(x))) > 0]
  }
  
} # end i

dimnames(incl.stages) <- list(cu_list$Conservation.Unit, stages)

#-----------------------------------------------------------------------------
# Write to Rds for call in calc-freshwater-exposure.r and summarize-freshwater-exposure.R
#-----------------------------------------------------------------------------

saveRDS(incl.stages, "output/freshwater-grid_included-cells.rds")
