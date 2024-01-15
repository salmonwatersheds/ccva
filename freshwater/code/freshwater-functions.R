###############################################################################
# Functions to load and summarize data relevant to exposure indicators
# 
# Contact: Stephanie Peacock <speacock@psf.ca>
# Date: April 17, 2023
###############################################################################

###############################################################################
# Function to load PCIC model output for given model and variable
###############################################################################

loadPCIC <- function(
    variable = "waterTemp", # Which variable to load? One of waterTemp or discharge
    model = "CanESM2", # Which GCM?
    rcp = 45 # Which emissions scenario? Options 45 or 85
){
  
  # Define location of data (external drive)
  root_dat <- "/Volumes/ClimateData/PCIC/hydro_model_out/fraser/"
  
  # Two variables in VICGL output
  varNames <- data.frame(
    fileName = c("waterTemp", "discharge"),
    Name = c("waterTemperature", "discharge")
  )
  
  v <- which(varNames$fileName == variable)
  
  #-----------------------------------------------------------------------------
  # There are different GCMs for the updated output
  # Want to load both discharge and waterTemp data for each model
  #-----------------------------------------------------------------------------
  
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
  z <- nc_open(paste0(root_dat, varNames$fileName[v], "_day_dynWat-VICGL_", gcms$modelName[j], "_rcp", rcp, "_", gcms$simulation[j], "_19450101-20991231_fraser.nc"))
  
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
  
  varExtracted0 <-  ncvar_get(z, varNames$Name[v])
  
  # Spatial variables
  lon <- ncvar_get(z, "lon")
  lat <- ncvar_get(z, "lat")
  
  # grid_points <- data.frame(
  #   id = c(1:(length(lon)*length(lat))),
  #   lon = rep(lon, length(lat)),
  #   lat = rep(lat, each = length(lon)),
  #   lon_id = rep(c(1:length(lon)), length(lat)),
  #   lat_id = rep(c(1:length(lat)), each = length(lon)))
  # write.csv(grid_points, file = "freshwater/output/PCIC-grid-points_bccoast.csv")
  
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
  spat.incl <- sort(read.csv("freshwater/output/PCIC_incl_overall.csv")$x)
  
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


#------------------------------------------------------------------------------
# Function to calculate 7-day average 
#------------------------------------------------------------------------------

calc_7DADM <- function(
    Tmax, # time-series of maximum daily temperatures
    date # associated dates
){
  
  if(sum(unique(diff(date)) != 1) > 0){
    stop("Not continuous time series")
    
  } else {
    
    n <- length(Tmax)
    
    dum <- rbind(
      c(rep(NA, 6), Tmax[1]),
      c(rep(NA, 5), Tmax[1:2]),
      c(rep(NA, 4), Tmax[1:3]),
      c(rep(NA, 3), Tmax[1:4]),
      c(rep(NA, 2), Tmax[1:5]),
      c(NA, Tmax[1:6]),
      cbind(
        Tmax[c(1:(n - 6))],
        Tmax[c(2:(n - 6 + 1))],
        Tmax[c(3:(n - 6 + 2))],
        Tmax[c(4:(n - 6 + 3))],
        Tmax[c(5:(n - 6 + 4))],
        Tmax[c(6:(n - 6 + 5))],
        Tmax[c(7:n)])
    )
    
    T7 <- apply(dum, 1, mean, na.rm = TRUE)
    return(T7)
  }
}

#------------------------------------------------------------------------------
# Function to calculate rolling 7-day maximum
#------------------------------------------------------------------------------
weeklyStat <- function(x, stat){
  n <- length(x)
  X <- rbind(
    x, 
    c(NA, x[1:(n-1)]),
    c(rep(NA, 2), x[c(1:(n - 2))]),
    c(rep(NA, 3), x[c(1:(n - 3))]),
    c(rep(NA, 4), x[c(1:(n - 4))]),
    c(rep(NA, 5), x[c(1:(n - 5))]),
    c(rep(NA, 6), x[c(1:(n - 6))])
  )    
  if(stat == "max"){
    xx <- apply(X, 2, max, na.rm = TRUE)
  } else if(stat == "mean"){
    xx <- apply(X, 2, mean, na.rm = TRUE)
  }

  return(xx)
    
}

#------------------------------------------------------------------------------
# Function to sum degree days outside of optimal range
#------------------------------------------------------------------------------
degdays <- function(
    x, # time series of temperature or flow
    thresh # vector of lower, upper thresholds
){
  
  sum(
    pmax(0, x - thresh[2], na.rm = TRUE), # Values above upper threshold
    pmax(0, thresh[1] - x, na.rm = TRUE))/30 # Values below lower threshold
  # divide by 30 years in period to get average annual dd
}

#------------------------------------------------------------------------------
# Function to Calculate days outside thresholds (d.o.t.)
#------------------------------------------------------------------------------
calcDOT <- function(x, # matrix of temperatures wiht dimension c(grid cells, days) 
                    thresholds # matrix of min/max temps for each grid cell
){
  dot <- numeric(dim(x)[1])  
  for(i in 1:(dim(x)[1])){
    dot[i] <- which(x[i, ] < thresholds[1, i] | x[i, ] > thresholds[2, i]) %>% length
  }
  return(dot)
}

#------------------------------------------------------------------------------
# Function to Calculate days below MAD thresholds 
#------------------------------------------------------------------------------
calcDBT <- function(x, # matrix of temperatures wiht dimension c(grid cells, days) 
                    mad_threshold # flow threshold (usually a % of mean annual discharge)
){
  dbt <- numeric(dim(x)[1])  
  for(i in 1:(dim(x)[1])){
    dbt[i] <- which(x[i, ] < mad_threshold[i]) %>% length
  }
  return(dbt)
}

#------------------------------------------------------------------------------
# Plot map
#------------------------------------------------------------------------------

map.stage <- function(
    cu_boundary.i,
    zoi.i,
    mig_paths.i,
    grid_polys = NA
){
  
  if(!is.na(zoi.i)){
    bounds <- cbind(st_bbox(zoi.i), st_bbox(mig_paths.i))
  bounds <- c(
    xmin = min(bounds[1, ]), 
    xmax = max(bounds[3,]), 
    ymin = min(bounds[2, ]), 
    ymax = max(bounds[4, ]))
  } else {
    bounds <- st_bbox(mig_paths.i)[c(1,3,2,4)]
  }
  
  par(bg = NA)
  
  plot(st_geometry(BC), border = NA, col = NA, axes = FALSE, las = 1, ylim = bounds[3:4], xlim = bounds[1:2], bty = "o")
  # plot(st_geometry(mig_paths.i), lwd = 3, col = cols[3], axes=TRUE, las = 1, ylim = bounds[3:4], xlim = bounds[1:2])
  # mtext(side = 3, line = 1, cus_to_keep$culabel[i])
  
  plot(grid_polys, add = TRUE, border = grey(0.8), col = NA)
  plot(BC, add = TRUE, col = NA, border = 1)
  
  
  if(!is.na(zoi.i)){
    plot(st_geometry(zoi.i), border = cols[1], col = paste0(cols[1], 50), lwd = 1.5, add = TRUE)
  }
  
  if(!is.na(cu_boundary.i)){
  plot(st_geometry(cu_boundary.i), border = 1, col = NA, lwd = 1, add = TRUE)
  }
  
  # coarser lakes and rivers
  plot(st_geometry(lakes0), border = 1, col = NA, add = TRUE)
  plot(st_geometry(rivers0), col = 1, add = TRUE)
  
  # legend("topleft", fill = c(paste0(cols[1], 50), cols[3], "white"), border = c(cols[1], cols[3], "#000000"), legend = c("Spawn ZOI", "Mig.", "CU boundary"), bg = "white")
  
}