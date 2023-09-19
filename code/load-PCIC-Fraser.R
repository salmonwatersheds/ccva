###############################################################################
# Function to load PCIC Updated model output
# 
# Contact: Stephanie Peacock <speacock@psf.ca>
# Date: Aug 31, 2023
###############################################################################
library(ncdf4)
library(dplyr)

###############################################################################
# Function to load PCIC model output for given model and variable
###############################################################################

loadPCIC <- function(
    variable = "waterTemp", # Which variable to load? One of waterTemp or discharge
    model = "CanESM2" # Which GCM?
    ){
  
  # Define location of data
  root_dat <- "data/raw-data/PCIC/hydro_model_out/fraser/"
  
  #-----------------------------------------------------------------------------
  # There are different GCMs for the updated output
  # Want to load both discharge and waterTemp data for each model
  #-----------------------------------------------------------------------------
  
  # Two variables in VICGL output
  varNames <- data.frame(
    fileName = c("waterTemp", "discharge"),
    Name = c("waterTemperature", "discharge")
  )
  
  i <- which(varNames$fileName == variable)
  
  # Vector of global climate model names (to extract netcdf output)
  gcms <- data.frame(
    modelName = c(
      "ACCESS1-0",
      "CanESM2",
      "CCSM4",
      "CNRM-CM5",
      "HadGEM2-ES",
      "MPI-ESM-LR"),
    simulation = c(
      "r1i1p1",
      "r1i1p1",
      "r2i1p1",
      "r1i1p1",
      "r1i1p1",
      "r3i1p1"
    )
  )
  j <- which(gcms$modelName == model)
  
  # Open netcdf file
  z <- nc_open(paste0(root_dat, varNames$fileName[i], "_day_dynWat-VICGL_", gcms$modelName[j], "_rcp45_", gcms$simulation[j], "_19450101-20991231_fraser.nc"))
  
  # # Explore dataset. What are the variables?
  # print(paste("File",z$filename,"contains",z$nvars,"variables"))
  # for( i in 1:z$nvars ) {
  #   v <- z$var[[i]]
  #   print(paste("Here is information on variable number",i))
  #   print(paste("   Name: ",v$name))
  #   print(paste("   Units:",v$units))
  #   print(paste("   Missing value:",v$missval))
  #   print(paste("   # dimensions :",v$ndims))
  #   print(paste("   Variable size:",v$varsize))
  # }
  
  # # Extract metadata
  # sink(paste0(root_dat, varNames$fileName[i], "_day_dynWat-VICGL_", gcms$modelName[j], "_rcp45_", gcms$simulation[j], "_19450101-20991231_fraser_metadata.txt"))
  # print(z)
  # sink()
  
  varExtracted0 <-  ncvar_get(z, varNames$Name[i])
  
  # Spatial variables
  lon <- ncvar_get(z, "lon")
  lat <- ncvar_get(z, "lat")

  # grid_points <- data.frame(
  #   id = c(1:(length(lon)*length(lat))),
  #   lon = rep(lon, length(lat)), 
  #   lat = rep(lat, each = length(lon)),
  #   lon_id = rep(c(1:length(lon)), length(lat)),
  #   lat_id = rep(c(1:length(lat)), each = length(lon)))
  # write.csv(grid_points, "output/PCIC-grid-points_Fraser.csv")
  
# Create time variable: all output runs from 19450101 to 20991231
  date <- as.Date(c(0:56612), origin = "1945-01-01")
  
  # Reduce dimensionality so that grid point[i, ] corresponds to Ts[i, ]
  varExtracted1 <- varExtracted0
  dim(varExtracted1) <- c(length(lon)*length(lat), length(date))
  
  #------------------------------------------------------------------------------
  # Dates since 1970-01-01 minus 7 days (for averaging)
  #------------------------------------------------------------------------------
  
  # Historical period is 1970-2000, so no need to keep data pre-1970
  date_keep <- which(date >= as.Date("1970-01-01") - 7)
  
  #------------------------------------------------------------------------------
  # Spatial grid cells used in analysis
  #------------------------------------------------------------------------------
  
  # PCIC included grid cells for each CU and life stage
  incl.stages <- readRDS("output/freshwater-grid_included-cells.rds")
  
  # Complete list of grid cells used in analysis
  spat.incl <- c()
  for(i in 1:n.CUs){
    for(j in 1:4){
      spat.incl <- sort(unique(c(spat.incl, incl.stages[[i,j]])))
    }
  } 
  
  #------------------------------------------------------------------------------
  # Trim to spatial and temporal extent needed
  #------------------------------------------------------------------------------
  varExtracted2 <- varExtracted1[spat.incl, date_keep] # approx 10% of the size

  dimnames(varExtracted2) <- list(
    spat.incl,
    as.character(date[date_keep])
  )
  
  rm(varExtracted0, varExtracted1)
  return(varExtracted2)
  
}

