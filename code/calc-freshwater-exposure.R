###############################################################################
# Code to summarize gridded hydrologic model output from the Pacific Climate 
# Impacts Consortium and calculate projected changes in exposure to stream
# temperature and flow outside of the historical norms for the given species 
# and life stage.
# 
# Contact: Stephanie Peacock <speacock@psf.ca>
# Date: Sept 13, 2023
###############################################################################

library(colorRamps)
library(mapproj)
library(proj4)
library(sf)
library(dplyr)
library(wesanderson)
library(ncdf4) # package for netcdf manipulation
library(zoo) # package with rollmean function

cols <- wes_palette("Darjeeling1")

source("code/functions.R")
source("code/load-PCIC-Fraser.R")

###############################################################################
# Load background spatial layers
###############################################################################

# Set broad bounds for clipping broad spatial layers
bounds <- c(xmin = -124.5, ymin = 50.85, xmax = -123.6, ymax = 51.7)

#------------------------------------------------------------------------------
# Lakes, rivers, and shorelines
#------------------------------------------------------------------------------

lakes <- readRDS("data/spatial/layers/waterbodies_250.rds")# %>% st_crop(bounds)
rivers <- readRDS("data/spatial/layers/watercourse_250.rds") #%>% st_crop(bounds)

bounds0 <- c(xmin = -127, ymin = 49, xmax = -116.5, ymax = 56)
lakes0 <- readRDS("data/spatial/layers/waterbodies_lowRes.rds")# %>% st_crop(bounds)
rivers0 <- readRDS("data/spatial/layers/watercourse_lowRes.rds") #%>% st_crop(bounds)
BC <- readRDS("data/spatial/layers/BC_lowRes.rds")


###############################################################################
# Define space and time variables for Fraser basin
###############################################################################

# Loading the PCIC data takes a minute or two
source("code/load-PCIC-prelimFraser.R")

# Variables created:
#  Ts = stream temperature (cell x date = 11794 x 56613)
#  grid_points = spatial data set of grid point centres
#  grid_polys = spatial polygons of each grid cell
#  date = object of class date with the days for each Ts grid cell

# Create time variable
# Note: PCIC output runs from 19450101 to 20991231, but is subsetting to date >= 1970-01-01 which is the earliest date used in the CCVAs (start of historical period)
date <- as.Date(0:47481, origin = "1970-01-01")

DOY <- as.numeric(strftime(date, format = "%j"))
months <- as.numeric(strftime(date, format = "%m"))
years <- as.numeric(strftime(date, format = "%Y"))

#-----------------------------------------------------------------------------
# Define periods for averages
#-----------------------------------------------------------------------------

periods <- data.frame(
  period = c("hist", "early", "mid", "late"),
  start = c(as.Date("1970-01-01"), as.Date("2010-01-01"), as.Date("2040-01-01"), as.Date("2070-01-01")),
  end = c(as.Date("1999-12-31"),as.Date("2039-12-31"),as.Date("2069-12-31"),as.Date("2099-12-31"))
)

periods$start_num <- match(periods$start, date)
periods$end_num <- match(periods$end, date)

#-----------------------------------------------------------------------------
# Define freshwater stages
#-----------------------------------------------------------------------------

stages <- c("adult_migration", "spawning", "eggs_alevin", "fw_rearing")

###############################################################################
# Population sensitivity data
###############################################################################

#------------------------------------------------------------------------------
# List of CUs in the PSE
#------------------------------------------------------------------------------

cu_list <- read.csv("data/cu_list.csv") %>% 
  subset(Species == "Lake sockeye" & Region == "Fraser" & Notes != "Extinct")

n.CUs <- length(unique(cu_list$Conservation.Unit))

#------------------------------------------------------------------------------
# Timing: Spawning
#------------------------------------------------------------------------------

# # Load other timing data from Sam (preliminary)
# timing <- read.csv("data/timing/Life_cycle_tables_20230731.csv") %>% subset(species == "SEL" & region == "fraser")
# timing$cuid <- cu_list$CUID[match(timing$culabel, cu_list$Conservation.Unit)]

timing0 <- read.csv("data/timing/4Life_cycle_table_by_CU.csv") %>% subset(species == "SEL" & region == "fraser")
timing0$cuid <- cu_list$CUID[match(timing0$culabel, cu_list$Conservation.Unit)]
timing0$culabel[is.na(timing0$cuid)] # CUs missing CUID are extinct
timing0 <- timing0[!is.na(timing0$cuid), ]

timing <- data.frame(
  culabel = rep(timing0$culabel, each = 4),
  cuid = rep(timing0$cuid, each = 4),
  stage = rep(c("eggs_alevin", "fw_rearing", "adult_migration", "spawning"), n.CUs),
  start = NA,
  end = NA
)

# Add timing start/end data
timing$start[timing$stage == "eggs_alevin"] <- round(timing0$inc_start)
timing$end[timing$stage == "eggs_alevin"] <- round(timing0$inc_end)
timing$start[timing$stage == "fw_rearing"] <- round(timing0$rear_start)
timing$end[timing$stage == "fw_rearing"] <- round(timing0$rear_end)
timing$start[timing$stage == "adult_migration"] <- round(timing0$mig_start)
timing$end[timing$stage == "adult_migration"] <- round(timing0$mig_end)
timing$start[timing$stage == "spawning"] <- round(timing0$sp_start)
timing$end[timing$stage == "spawning"] <- round(timing0$sp_end)

write.csv(timing, "output/freshwater_timing_FraserSEL.csv", row.names = FALSE)
# # Adjust stage labels to match
# timing$stage[timing$stage == "Incubation"] <- "eggs_alevin"
# timing$stage[timing$stage == "Rearing"] <- "fw_rearing"
# timing$stage[timing$stage == "Upriver Migration"] <- "adult_migration"
# timing$stage[timing$stage == "Spawning"] <- "spawning"
# 


# # **Initial focus on Fraser CUs with specific thermal tolerance data **
# cus_to_keep <-data.frame(
#   culabel = c("Chilko-Early Summer", "Chilko-Summer", "Takla-Trembleur-Early Stuart (cyclic)",  "Francois-Fraser-Summer", "Nadina-Francois-Early Summer", "Quesnel-Summer (cyclic)",
#                  "Anderson-Seton-Early Summer",  "Kamloops-Early Summer", "Shuswap-Early Summer (cyclic)", "Shuswap-Late (cyclic)",  "Harrison-Upstream Migrating-Late", "Harrison-Downstream Migrating-Late"),
#   population = c(rep("Chilko", 2), "EarlyStuart", rep("Nechako", 2), "Quesnel", "Gates", rep("Lower Adams", 3), "Weaver", "Harrison"))
# cus_to_keep$cuid <- cu_list$CUID[match(cus_to_keep$culabel, cu_list$Conservation.Unit)]
# popCol <- c(Chilko = "#384E9D", Nechako = "#EF1F25", Gates = "#EC7E21", EarlyStuart = "#7C2781", Quesnel = "#0B803F", Weaver = "#6ACADC", Harrison = "#FF40FF")

#------------------------------------------------------------------------------
# Spawning Zone of Influence: Fraser SEL
#------------------------------------------------------------------------------

spawn_zoi <- st_read(dsn = "data/spatial/ZOI/fraser-spawning-zoi/fraser_spawning_zoi_SEL.shp") %>% st_transform(crs = 4269)

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


###############################################################################
# Global climate model names
###############################################################################
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

n.models <- length(gcms$modelName)

###############################################################################
# Loop through models, life stages, CUs
###############################################################################

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------

fw_output <- array(NA,
                   dim = c(n.models, n.CUs, length(stages), 4, 4, 3),
                   dimnames = list(
                     gcms$modelName,
                     cu_list$Conservation.Unit, 
                     stages,
                     c("optimalTemp", "criticalTemp", "optimalFlow", "criticalFlow"),
                     c("hist", "early", "mid", "late"), 
                     c("median", "lower", "upper")))

# Create list to store full spatial output for plotting maps
# Later will calculate median across all GCMs
fw_spat <- list(); length(fw_spat) <- n.CUs*length(stages)
dim(fw_spat) <- c( n.CUs, length(stages))
dimnames(fw_spat) <- list(cu_list$Conservation.Unit, stages)


# Store temperature thresholds for each stage
tempThresh <- data.frame(
  cuid = rep(cu_list$CUID, each = length(stages)*n.models),
  culabel = rep(cu_list$Conservation.Unit, each = length(stages)*n.models),
  stage = rep(rep(stages, each = n.models), n.CUs),
  model = rep(gcms$modelName, length(stages)*n.CUs),
  critical_max = rep(NA, n.models*n.CUs*length(stages)),
  optimal_max = rep(NA, n.models*n.CUs*length(stages)))

# Track time to run (takes ~ 20 mins for Fraser SEL CUs)
start.time <- Sys.time()
lap.time <- Sys.time()
time.for.models <- matrix(NA, nrow = n.CUs, ncol = n.models)

#------------------------------------------------------------------------------
# Loop through all Global Climate Models
#------------------------------------------------------------------------------

for(m in 1:n.models){
  
  # Load model output; this takes a minute
  Ts <- loadPCIC(
    variable = "waterTemp", # Which variable to load? One of waterTemp or discharge
    model = gcms$modelName[m] # Which GCM?
  )
  
  Qs <- loadPCIC(
    variable = "discharge",
    model = gcms$modelName[m]
  )
  
  # Create grid points; same for all models so only need to do once
  if(m == 1){
    lon <- as.numeric(dimnames(Ts)[[1]])
    lat <- as.numeric(dimnames(Ts)[[2]])
    
    grid_points <- data.frame(
      id = c(1:(length(lon)*length(lat))),
      lon = rep(lon, length(lat)), 
      lat = rep(lat, each = length(lon))) %>%
      st_as_sf(coords = c("lon", "lat"), crs = 4269)
    
    d <- 1/16 # spacing of grid points
  }
  
  # Reduce dimensionality so that grid point[i, ] corresponds to Ts[i, ]
  dim(Ts) <- c(length(lon)*length(lat), length(date))
  dim(Qs) <- c(length(lon)*length(lat), length(date))
  
  # Change Ts from degrees Kelvin to Celsius
  if(min(Ts[, 1000], na.rm = TRUE) > 200){
    Ts <- Ts - 273.15
  }
  #------------------------------------------------------------------------------
  # Loop through each CU
  #------------------------------------------------------------------------------
  
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
    
    #------------------------------------------------------------------------------
    # Loop through each life stage
    #------------------------------------------------------------------------------
    
    for(j in 1:length(stages)){
      
      # If there are timing data
      if(length(which(timing$cuid == cuid & timing$stage == stages[j])) == 0){
        warning(paste0("No timing data for ", cu_list$Conservation.Unit[i], " ", stages[j]))
        
      } else {
        
        # Extract timing and covert to day-of-year
        timing.ij_prelim1 <- timing[which(timing$cuid == cuid & timing$stage == stages[j])[1], c("start", "end")] %>% 
          as.numeric()
        if(timing.ij_prelim1[2] - timing.ij_prelim1[1] + 1 == 365){
          timing.ij <- c(1, 365)
        } 
        
        #-----------------------------------------------------------------------------
        # Identify grid cells for given life stage
        #-----------------------------------------------------------------------------
        # For adult migration and freshwater rearing, include downstream migration route
        if(stages[j] %in% c("adult_migration", "fw_rearing")){
          intrscts <- st_is_within_distance(grid_points, mig_paths.i, 3000, sparse = FALSE)  
          incl <- which(c(intrscts[, 1]) == TRUE)
          # There may be more than one line segment
          if(dim(intrscts)[2] > 1){
            for(k in 2:dim(intrscts)[2]){
              incl <- unique(c(incl, which(c(intrscts[, k]) == TRUE)))
            }
          }
        } 
        
        # For spawning, eggs_alevin, use spawning ZOI
        if(stages[j] %in% c("spawning", "eggs_alevin")){
          intrscts <- st_intersects(grid_points, zoi.i, sparse = FALSE)
          incl <- which(c(intrscts) == TRUE)
        }
        
        # For freshwater rearing, include rearing lake (CU boundary)
        if(stages[j] == "fw_rearing"){
          
          # CU boundary is a multipolygon and was throwing error: 
          p <- st_make_valid(cu_boundary.i)
          
          intrscts <- st_intersects(grid_points, p, sparse = FALSE)
          incl <- unique(c(incl, which(c(intrscts) == TRUE)))
        } 
        
        # Exclude cells that have no data
        incl <- incl[apply(Ts[incl, ], 1, function(x) sum(!is.na(x))) > 0]
        
        # Highlight those cells
        
        grid_polys_incl <- st_as_sf(data.frame(
          id = rep(grid_points$id[incl], each = 5),
          rep(c("SW0", "NW", "NE", "SE", "SW1"), length(incl)),
          lon = c(rep(rep(lon, length(lat))[incl], each = 5) + rep(c(- d/2, -d/2, d/2,  d/2, -d/2), length(incl))),
          lat = c(rep(rep(lat, each = length(lon))[incl], each = 5) + rep(c(- d/2,  d/2, d/2, -d/2, -d/2), length(incl)))
        ), coords = c("lon", "lat"), crs = 4269) %>% 
          group_by(id) %>%
          summarise(geometry = st_combine(geometry)) %>%
          st_cast("POLYGON") 
        
        # Dimensions
        n.grid <- length(incl)
        
        # Setup array to hold full spatial output
        fw_spat[[i, j]] <- array(data = NA,
                                 dim = c(n.models, 4, 4, n.grid),
                                 dimnames = list(
                                   gcms$modelName,
                                   c("optimalTemp", "criticalTemp", "optimalFlow", "criticalFlow"),
                                   c("hist", "early", "mid", "late"),
                                   incl)
        )
        
        #-----------------------------------------------------------------------------
        # Calculate rolling weekly max (Ts) or mean (Qs) over life-stage period
        #-----------------------------------------------------------------------------
        
        # Select dates in spawn timing, extending earlier 7 days for calculating 7-day rolling avg
        periods.ind <- rbind(
          hist = which(DOY %in% c(c(timing.ij[1] - 7):timing.ij[2]) & years %in% c(1970:1999)),
          early = which(DOY %in% c(c(timing.ij[1] - 7):timing.ij[2]) & years %in% c(2010:2039)),
          mid = which(DOY %in% c(c(timing.ij[1] - 7):timing.ij[2]) & years %in% c(2040:2069)),
          late = which(DOY %in% c(c(timing.ij[1] - 7):timing.ij[2]) & years %in% c(2070:2099)))
        
        # Dimensions
        n.days <- ncol(periods.ind) - 7
        
        
        # Create arrays
        Ts.weeklyMax <- array(
          NA, 
          dim = c(4, n.grid, n.days), 
          dimnames = list(c("hist", "early", "mid", "late")))
        Qs.weeklyMean <- Ts.weeklyMax
        
        # Loop through each period
        for(p in 1:4){ 
          Ts.weeklyMax[p, , ] <- t(apply(Ts[incl, periods.ind[p,]], 1, weeklyStat, stat = "max"))[, 8:ncol(periods.ind)]
          Qs.weeklyMean[p, , ] <- t(apply(Qs[incl, periods.ind[p,]], 1, weeklyStat, stat = "mean"))[, 8:ncol(periods.ind)]
        }
        
        # Checks in case -Inf (not getting this any more?)
        Ts.weeklyMax[which(Ts.weeklyMax == -Inf)] <- NA
        Qs.weeklyMean[which(Qs.weeklyMean == -Inf)] <- NA
        
        
        #-----------------------------------------------------------------------------
        # Stream temperature 
        #-----------------------------------------------------------------------------
        
        # Historical temperature
        #  # dimensions: space (grid cells) x time (days in 30-year period)
        
        # Calculate acute (critical) and chronic (optimal) temperature thresholds
        thresh <- apply(Ts.weeklyMax[1, , ], 1, quantile, c(0.1, 0.9, 0.025, 0.975), na.rm = TRUE)
        rownames(thresh) <- c("opt_min", "opt_max", "crit_min", "crit_max")
        
        if(3 == 2){
          par(mfrow = c(1,1), mar = c(4,4,2,1), oma = c(0,0,2,0))
          hist(Ts.weeklyMax[1, 10, ], las = 1, xlab = "Weekly max. (˚C)", main = "", col = paste0(cols[5], 60), border = NA)
          abline(v = thresh[c("opt_min", "opt_max"), 10], col = cols[2], lty = 2)
          abline(v = thresh[c("crit_min", "crit_max"), 10], col = cols[2])
          for(k in 2:4){
            
            hist(Ts.weeklyMax[k, 10, ], las = 1, xlab = "", main = "", col = paste0(cols[c(3,4,1)[k-1]], 60), border = NA, add = TRUE)
          }
          mtext(side = 3, outer = TRUE, paste(cu_list$culabel[i], stages[j]))
          legend("topright", fill = paste0(cols[c(5,3,4,1)], 60), border = NA, legend = c("hist", "early", "mid", "late"))
        }
        
        # Apply hard thresholds where historical temps are extreme
        # ** Might need to think about also imposing higher crit_max for northern CUs where historical temp distribution is cold...
        thresh["crit_max", which(thresh["crit_max", ] > 26)] <- 26 # Critical max 26˚C, Richter & Kolmes (2005)
        thresh["opt_max", which(thresh["opt_max", ] > 24)] <- 24 # Chronic max 24˚C, Sullivan 2000
        
        #  ** Ignore minimum thresholds for now **
        # thresh["opt_min", which(thresh["opt_min", ] < c(8, 10, 8, 7)[j])] <- c(8, 10, 8, 7)[j] # Example of minimum temps from Mayer (in review)
        thresh["opt_min", ] <- 0
        thresh["crit_min", ] <- 0
        
        # Check that optimal is always inside critical
        if(length(which(thresh["crit_max", ] < thresh["opt_max", ])) > 0){
          stop("crit_max < opt_max")
        }
        # if(length(which(thresh["crit_min", ] > thresh["opt_min", ])) > 0){
        #   stop("crit_min > opt_min")
        # }
        
        # Store median thresholds for output
        tempThresh[which(tempThresh$cuid == cuid & tempThresh$stage == stages[j]), c("critical_min", "critical_max")] <- apply(thresh[c("crit_min", "crit_max"), ], 1, median)
        tempThresh[which(tempThresh$cuid == cuid & tempThresh$stage == stages[j]), c("optimal_min", "optimal_max")] <- apply(thresh[c("opt_min", "opt_max"), ], 1, median)
        
        # Calculate days outside of threshold
        Ts.dot <- array(NA, 
                        dim = c(2, 4, length(incl)),
                        dimnames = list(c("opt", "crit"), c("hist", "early", "mid", 'late'), NULL)
        )
        for(p in 1:4){
          for(k in 1:2){
            Ts.dot[k, p, ] <- calcDOT(x = Ts.weeklyMax[p, , ], 
                                      thresholds = thresh[list(c("opt_min", "opt_max"), c("crit_min", "crit_max"))[[k]], ])
          }}
        
        # Calculate % of days outside window
        Ts.perc <- Ts.dot/n.days
        
        # Store full spatial output
        fw_spat[[i,j]][m, c("optimalTemp", "criticalTemp"), , ] <- Ts.perc
        
        # store summary output
        for(p in 1:4){ # for each period
          for(k in 1:2){ # for optimal and critical
            # Store output in end array
            fw_output[m, i, j, c("optimalTemp", "criticalTemp")[k], p, ] <- c(
              median(Ts.perc[k, p, ]),
              range(Ts.perc[k, p, ]))
            
          }}
        
        # # Calculate % change in the number of days outside window relative to historical
        # Ts.pchange <- array(NA, 
        #                     dim = c(2, 3, 3), 
        #                     dimnames = list(c("acute", "chronic"), 
        #                                     c("early", "mid", "late"), 
        #                                     c("median", "min", "max"))
        # )
        # 
        # for(p in 1:3){
        #   for(k in 1:2){
        #     perc <- (Ts.dot[k, p + 1, ] - Ts.dot[k, 1, ])/Ts.dot[k, 1, ] * 100
        #     Ts.pchange[k, p, "median"] <- median(perc)
        #     Ts.pchange[k, p, c("min", "max")] <- range(perc)
        #     
        #     # Store output in end array
        #     fw_exposure[k, i, j, p, ] <- Ts.pchange[k, p, ]
        # }}
        
        
        #-----------------------------------------------------------------------------
        # Flow
        #-----------------------------------------------------------------------------
        
        # Calculate MAD over historical period for each included grid cell
        mad <- apply(Qs[incl, c(periods[periods$period == 'hist', 'start_num']:periods[periods$period == "hist", "end_num"])], 1, mean, na.rm = TRUE)
        
        # Calculate acute (critical) and chronic (optimal) temperature thresholds
        thresh_mad <- mad %*% matrix(c(0.2, 0.1), nrow = 1)
        colnames(thresh_mad) <- c("opt", "crit")
        
        # Calculate days below % mad threshold
        Qs.dbt <- array(NA, 
                        dim = c(2, 4, n.grid),
                        dimnames = list(c("opt", "crit"), 
                                        c("hist", "early", "mid", 'late'), 
                                        NULL)
        )
        for(p in 1:4){
          for(k in 1:2){
            Qs.dbt[k, p, ] <- calcDBT(x = Qs.weeklyMean[p, , ], 
                                      mad_threshold = c(thresh_mad[, c("opt", "crit")[k]]))
          }}
        
        
        # Calculate % of days outside window
        Qs.perc <- Qs.dbt/n.days
        
        # Store full spatial output
        fw_spat[[i,j]][m, c("optimalFlow", "criticalFlow"), , ] <- Qs.perc
        
        for(p in 1:4){ # for each period
          for(k in 1:2){
            # Store output in end array
            fw_output[m, i, j, c("optimalFlow", "criticalFlow")[k], p, ] <- c(
              median(Qs.perc[c("opt", "crit")[k], p, ]),
              range(Qs.perc[c("opt", "crit")[k], p, ]))
          }}
        
      } # end if timing data
      print(".")
      
    } # end life stage j
    
    print(paste0("Done ", cu_list$Conservation.Unit[i]))
    
  } # end CU i
  
} # end model m
print(Sys.time() - start.time)
# Time difference of 7.366207 hours
###############################################################################
# Save output
###############################################################################

write.csv(tempThresh, file = "output/tempThresholds_FraserSEL.csv")
saveRDS(fw_spat, file = "output/freshwater_spat_FraserSEL.rds")
saveRDS(fw_output, file = "output/freshwater_output_FraserSEL.rds")

write.csv(grid_points, file = "output/PCIC_grid-points_Fraser.csv")
