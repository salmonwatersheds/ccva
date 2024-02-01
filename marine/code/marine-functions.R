###############################################################################
# Functions to load and summarize data relevant to exposure indicators
# 
# Contact: Stephanie Peacock <speacock@psf.ca>
# Date: April 17, 2023
###############################################################################
library(abind)

###############################################################################
# Function to load NOAA model output for given model and variable
###############################################################################

loadNOAA <- function(
    variable = "SST", # Which variable to load? One of "PP", "pH", "SSS", "SST"
    model = "Can-ESM2", # Which GCM? One of "ACCESS1-0", "Can-ESM2", "CCSM4", "CNRM-CM5", "HADGEM2-ES", "MPI-ESM-LR"
    rcp = 45, # Which emissions scenario? Options 45 or 85
    spat_id = NULL # ids of spatial grid cells to keep ("marine/output/incl_NOAA_EM.csv" or "marine/output/incl_NOAA_marRear.csv")
    ){
  
  # Define location of data (external drive)
  root_dat <- "/Volumes/ClimateData/NOAA/"
  
  #-----------------------------------------------------------------------------
  # There are different GCMs for the updated output
  #-----------------------------------------------------------------------------
  
  # Variables in NOAA data dump
  #   tos is the temperature at surface
  #   sos is the sea surface salinity
  #   ph is the surface ph
  #   intpp is the vertically integrated net primary carbon productivity rate from phytoplankton
  
  varNames <- data.frame(
    fileName = c("intpp", "ph", "sos", "tos"),
    Name = c("PP", "pH", "SSS", "SST")
  )
  
  i <- which(varNames$Name == variable)
  
  # # Vector of global climate model names (to extract netcdf output)
  # gcms <- data.frame(
  #   modelName = c(
  #     "ACCESS1-0",
  #     "Can-ESM2",
  #     "CCSM4",
  #     "CNRM-CM5",
  #     "HADGEM2-ES",
  #     "MPI-ESM-LR"),
  # 
  # )
  
  # Open netcdf file
  z <- nc_open(paste0(root_dat, varNames$fileName[i], "_Omon_", model, "_rcp", rcp, "_r1i1p1_200601-210012.1x1.nc"))
  
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
  
  # Get projected
  varExtracted_proj <-  ncvar_get(z, varNames$fileName[i])
  
  # Get historical
  z_hist <- nc_open(paste0(root_dat, varNames$fileName[i], "_Omon_", model, "_historical_r1i1p1_190101-200512.1x1.nc"))
  varExtracted_hist <-  ncvar_get(z_hist, varNames$fileName[i])
  
  varExtracted <- abind(varExtracted_hist, varExtracted_proj, along = 3)
  
  dates <- as.Date(paste(01, rep(1:12, 200), rep(1901:2100, each = 12), sep = "-"), format = "%d-%m-%Y")
  
  # Spatial variables
  lon <- ncvar_get(z, "lon")
  lat <- ncvar_get(z, "lat")

  # grid_points <- data.frame(
  #   id = c(1:(length(lon)*length(lat))),
  #   lon = rep(lon, length(lat)),
  #   lat = rep(lat, each = length(lon)),
  #   lon_id = rep(c(1:length(lon)), length(lat)),
  #   lat_id = rep(c(1:length(lat)), each = length(lon)))
  # 
  # grid_points$hasData <- rep(0, dim(grid_points)[1])
  # grid_points$hasData[notNA] <- 1
  #  write.csv(grid_points, "output/NOAA-grid-points.csv")

  # Reduce dimensionality so that grid point[i, ] corresponds to Ts[i, ]
  varExtracted1 <- varExtracted
  dim(varExtracted1) <- c(length(lon)*length(lat), length(dates))
  
   
  #------------------------------------------------------------------------------
  # Trim to relevant spatial extent; if not given then trim where there's no data
  #------------------------------------------------------------------------------
  
  if(!missing(spat_id)){
    spat_keep <- spat_id
  } else {
    spat_keep <- which(!is.na(varExtracted1[, 1]))
  }
  
  spat_keep <- spat_keep[which(!is.na(varExtracted1[spat_keep, 1]))]
  
  #------------------------------------------------------------------------------
  # Dates since 1970-01-01 minus 7 days (for averaging)
  #------------------------------------------------------------------------------
  
  # Historical period is 1970-2000, so no need to keep data pre-1970
  date_keep <- which(dates >= as.Date("1970-01-01"))
  
  
  #------------------------------------------------------------------------------
  # Trim to spatial and temporal extent needed
  #------------------------------------------------------------------------------
  varExtracted2 <- varExtracted1[spat_keep, date_keep] # approx 10% of the size

  dimnames(varExtracted2) <- list(
    spat_keep,
    as.character(dates[date_keep])
  )
  
  rm(varExtracted, varExtracted1)
  
  # Convert K to C
  if(min(varExtracted2[1,1:100], na.rm = TRUE) > 200 & i == 4){
    varExtracted2 <- varExtracted2 - 273.15
  }
  return(varExtracted2)
  
}

