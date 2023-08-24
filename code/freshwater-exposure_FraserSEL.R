###############################################################################
# Code to summarize gridded hydrologic model output from the Pacific Climate 
# Impacts Consortium and calculate projected changes in exposure to stream
# temperatures outside of the optimal range for the given species and life stage.
# 
# Contact: Stephanie Peacock <speacock@psf.ca>
# Date: April 17, 2023
###############################################################################

library(colorRamps)
library(mapproj)
library(proj4)
library(sf)
# library(leaflet)
library(dplyr)
library(wesanderson)
library(ncdf4) # package for netcdf manipulation
library(zoo) # package with rollmean function

cols <- wes_palette("Darjeeling1")

source("code/functions.R")
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
# PCIC *Preliminary* climate model output for the Fraser Basin
###############################################################################

# Loading the PCIC data takes a minute or two
source("code/load-PCIC-prelimFraser.R")

# Variables created:
#  Ts = stream temperature (cell x date = 11794 x 56613)
#  grid_points = spatial data set of grid point centres
#  grid_polys = spatial polygons of each grid cell
#  date = object of class date with the days for each Ts grid cell

DOY <- as.numeric(strftime(date, format = "%j"))
months <- as.numeric(strftime(date, format = "%m"))
years <- as.numeric(strftime(date, format = "%Y"))

#-----------------------------------------------------------------------------
# Define periods for averages
#-----------------------------------------------------------------------------

periods <- data.frame(
  period = c("hist", "early", "mid", "late"),
  start = c(as.Date("1970-01-01"), as.Date("2010-01-01"), as.Date("2040-01-01"), as.Date("2070-01-01")),
  end = c(as.Date("2000-01-01"),as.Date("2040-01-01"),as.Date("2070-01-01"),as.Date("2100-01-01"))
)

periods$start_num <- match(periods$start, date)
periods$end_num <- match(periods$end - 1, date)

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

cu_list <- read.csv("https://raw.githubusercontent.com/salmonwatersheds/tech-report/main/tables/appendix1.csv") %>% subset(Species == "Lake sockeye" & Region == "Fraser")

# #-----------------------------------------------------------------------------
# # Optimal temperature ranges
# # Environmental Protection Division. 2001. Water Quality Guidelines for Temperature: Overview Report. Available from https://www2.gov.bc.ca/assets/gov/environment/air-land-water/water/waterquality/water-quality-guidelines/approved-wqgs/temperature-or.pdf.
# #-----------------------------------------------------------------------------
# 
# topt.all <- read.csv("1-exposure/data/optimum-temps_BC2001.csv")
# topt <- as.numeric(topt.all[which(topt.all$Species == "Sockeye" & topt.all$Life.stage == "Migration"), c("min", "max")])
# topt2 <- read.csv("1-exposure/data/topt_Eliason.csv")


#------------------------------------------------------------------------------
# Timing: Spawning
#------------------------------------------------------------------------------

# timing <- read.csv("data/timing/spawntiming_byCU.csv")
# 
# timing$spawn_start <- timing$start_spawn_025
# timing$spawn_end <- timing$end_spawn_975
# 
# # Pull in PSE display name for CU
# timing$culabel <- cu_list$culabel[match(timing$cuid, cu_list$cuid)]

# Load other timing data from Sam (preliminary)
timing <- read.csv("data/timing/Life_cycle_tables_20230731.csv") %>% subset(species == "SEL" & region == "fraser")
timing$cuid <- cu_list$CUID[match(timing$culabel, cu_list$Conservation.Unit)]

# Adjust stage labels to match
timing$stage[timing$stage == "Incubation"] <- "eggs_alevin"
timing$stage[timing$stage == "Rearing"] <- "fw_rearing"
timing$stage[timing$stage == "Upriver Migration"] <- "adult_migration"
timing$stage[timing$stage == "Spawning"] <- "spawning"


# **Initial focus on Fraser CUs with specific thermal tolerance data **
cus_to_keep <-data.frame(
  culabel = c("Chilko-Early Summer", "Chilko-Summer", "Takla-Trembleur-Early Stuart (cyclic)",  "Francois-Fraser-Summer", "Nadina-Francois-Early Summer", "Quesnel-Summer (cyclic)",
                 "Anderson-Seton-Early Summer",  "Kamloops-Early Summer", "Shuswap-Early Summer (cyclic)", "Shuswap-Late (cyclic)",  "Harrison-Upstream Migrating-Late", "Harrison-Downstream Migrating-Late"),
  population = c(rep("Chilko", 2), "EarlyStuart", rep("Nechako", 2), "Quesnel", "Gates", rep("Lower Adams", 3), "Weaver", "Harrison"))
cus_to_keep$cuid <- cu_list$CUID[match(cus_to_keep$culabel, cu_list$Conservation.Unit)]
popCol <- c(Chilko = "#384E9D", Nechako = "#EF1F25", Gates = "#EC7E21", EarlyStuart = "#7C2781", Quesnel = "#0B803F", Weaver = "#6ACADC", Harrison = "#FF40FF")

#------------------------------------------------------------------------------
# Spawning Zone of Influence: Fraser SEL
#------------------------------------------------------------------------------

spawn_zoi <- st_read(dsn = "data/spatial/ZOI/fraser-spawning-zoi/fraser_spawning_zoi_SEL.shp") %>% st_transform(crs = 4269)

#------------------------------------------------------------------------------
# Adult migration lines: Fraser SEL
#------------------------------------------------------------------------------
mig_paths <- st_read(dsn = "data/spatial/fw-migration/spawn_timing_migration_paths_wCU_fraser.shp") %>% st_transform(crs = 4269)

# Check these paths
plot(mig_paths)

# #------------------------------------------------------------------------------
# # Spawning points and lines: Fraser SEL
# #------------------------------------------------------------------------------
# 
# spawning_points <- st_read(dsn = "1-exposure/data/spatial/spawning_points_fraser/spawning_points_fraser.shp") %>% st_transform(crs = 4269) %>% st_crop(bounds0)
# 
# spawning_lines <- st_read(dsn = "1-exposure/data/spatial/spawning_lines_fraser/spawning_lines_fraser.shp")%>% st_transform(crs = 4269) %>% st_crop(bounds0) 

#------------------------------------------------------------------------------
# Conservation Unit Boundaries: SEL
#------------------------------------------------------------------------------

shp_path <- "data/spatial/CU-boundaries/Lake_Type_Sockeye_Salmon_CU_Boundary/SEL_CU_BOUNDARY_En"
cu_boundary <- st_read(dsn = paste0(shp_path, ".shp"), promote_to_multi = FALSE)

###############################################################################
# Select CU for analysis
###############################################################################

# cuids <- timing$cuid[which(timing$region == "Fraser" & timing$SPECIES_QUALIFIED == "SEL" & !is.na(timing$cuid))] # 23 Fraser SEL Cus
cuids <- cus_to_keep$cuid

# Setup arrays to store output
fw_exposure <- array(NA,
            dim = c(4, length(cuids), length(stages), 4, 3),
            dimnames = list(
              c("temp_opt", "temp_crit", "lowflow_opt", "lowflow_crit"),
              cus_to_keep$culabel, 
              stages,
              c("hist", "early", "mid", "late"), 
              c("median", "lower", "upper")))

# Create list to store full spatial output for plotting maps
fw_exposure.spat <- list(); length(fw_exposure.spat) <- length(cuids)*length(stages) 
dim(fw_exposure.spat) <- c(length(cuids), length(stages))

# Store temperature thresholds for each stage
tempThresh <- data.frame(
  cuid = rep(cus_to_keep$cuid, each = length(stages)),
  culabel = rep(cus_to_keep$culabel, each = length(stages)),
  stage = rep(stages, length(cuids)),
  critical_min = rep(NA, length(cuids)*length(stages)),
  critical_max = rep(NA, length(cuids)*length(stages)),
  optimal_min = rep(NA, length(cuids)*length(stages)),
  optimal_max = rep(NA, length(cuids)*length(stages)))
#------------------------------------------------------------------------------
# For each CU
#------------------------------------------------------------------------------
# Track time to run (takes ~ 20 mins for Fraser SEL CUs)
start.time <- Sys.time()

# Loop through each CU
for(i in 1:length(cuids)){
  
  # Subset temporal and spatial data
  cu_boundary.i <- cu_boundary[which(cu_boundary$FULL_CU_IN == cu_list$Full.CU.Index[which(cu_list$CUID == cuids[i])]),]
  
  zoi.i <- spawn_zoi[which(spawn_zoi$cuid == cuids[i]), ]
  
  mig_paths.i <- mig_paths[which(mig_paths$cuid == cuids[i]), ]
  
  #------------------------------------------------------------------------------
  # For each life stage
  #------------------------------------------------------------------------------
  for(j in 1:length(stages)){
    
    # If there are timing data
    if(length(which(timing$cuid == cuids[i] & timing$stage == stages[j])) > 0){
      
      # Extract timing and covert to day-of-year
      timing.ij_prelim1 <- timing[which(timing$cuid == cuids[i] & timing$stage == stages[j])[1], c("start", "end")] %>% 
      as.character() %>%
      as.Date(format = "%Y-%m-%d")
      if(timing.ij_prelim1[2] - timing.ij_prelim1[1] + 1 == 365){
        timing.ij <- c(1, 365)
      } else {
        timing.ij <- strftime(timing.ij_prelim1, format = "%j") %>% 
          as.numeric() %>%
          round()
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
      
      # Highlight those cells
      grid_polys_incl <- st_as_sf(data.frame(
        V2 = rep(spat_grid$V2[incl], each = 5),
        rep(c("SW0", "NW", "NE", "SE", "SW1"), length(incl)),
        lon = c(rep(spat_grid$lon[incl], each = 5) + rep(c(- d/2, -d/2, d/2,  d/2, -d/2), length(incl))),
        lat = c(rep(spat_grid$lat[incl], each = 5) + rep(c(- d/2,  d/2, d/2, -d/2, -d/2), length(incl)))
      ), coords = c("lon", "lat"), crs = 4269) %>% 
        group_by(V2) %>%
        summarise(geometry = st_combine(geometry)) %>%
        st_cast("POLYGON") 
      
      # Dimensions
      n.grid <- length(incl)
      
      # Setup array to hold full spatial output
      fw_exposure.spat[[i,j]] <- array(data = NA,
                              dim = c(4, 4, n.grid),
                              dimnames = list(
                                c("temp_opt", "temp_crit", "lowflow_opt", "lowflow_crit"),
                                c("hist", "early", "mid", "late"),
                                NULL)
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
      
      hist(Ts.weeklyMax[k, , ], las = 1, xlab = "", main = "", col = paste0(cols[c(3,4,1)[k-1]], 60), border = NA, add = TRUE)
    }
    mtext(side = 3, outer = TRUE, paste(cus_to_keep$culabel[i], stages[j]))
    legend("topright", fill = paste0(cols[c(5,3,4,1)], 60), border = NA, legend = c("hist", "early", "mid", "late"))
  }
  
  # Apply hard thresholds where historical temps are extreme
  # ** Might need to think about also imposing higher crit_max for northern CUs where historical temp distribution is cold...
  thresh["crit_max", which(thresh["crit_max", ] > 26)] <- 26 # Critical max 26˚C, Richter & Kolmes (2005)
  thresh["opt_max", which(thresh["opt_max", ] > 24)] <- 24 # Chronic max 24˚C, Sullivan 2000
  
  thresh["opt_min", which(thresh["opt_min", ] < c(8, 10, 8, 7)[j])] <- c(8, 10, 8, 7)[j] # Example of minimum temps from Mayer (in review)
  
  # Check that optimal is always inside critical
  if(length(which(thresh["crit_max", ] < thresh["opt_max", ])) > 0){
    stop("crit_max < opt_max")
  }
  if(length(which(thresh["crit_min", ] > thresh["opt_min", ])) > 0){
    stop("crit_min > opt_min")
  }
  
  # Store median thresholds for output
  tempThresh[which(tempThresh$cuid == cuids[i] & tempThresh$stage == stages[j]), c("critical_min", "critical_max")] <- apply(thresh[c("crit_min", "crit_max"), ], 1, median)
  tempThresh[which(tempThresh$cuid == cuids[i] & tempThresh$stage == stages[j]), c("optimal_min", "optimal_max")] <- apply(thresh[c("opt_min", "opt_max"), ], 1, median)
  
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
  fw_exposure.spat[[i,j]][c("temp_opt", "temp_crit"), , ] <- Ts.perc
    
    # store summary output
  for(p in 1:4){ # for each period
      for(k in 1:2){
        # Store output in end array
        fw_exposure[k, i, j, p, ] <- c(
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
  
  #--------------
  # Calculate MAD over historical period for each included grid cell
  #-------------
  
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
  fw_exposure.spat[[i,j]][c("lowflow_opt", "lowflow_crit"), , ] <- Qs.perc
  
  for(p in 1:4){ # for each period
    for(k in 1:2){
      # Store output in end array
      fw_exposure[c("lowflow_opt", "lowflow_crit")[k], i, j, p, ] <- c(
        median(Qs.perc[c("opt", "crit")[k], p, ]),
        range(Qs.perc[c("opt", "crit")[k], p, ]))
    }}
  
    } # end if timing data
    print(".")
  } # end life stage j
  print(paste0("Done ", cus_to_keep$culabel[i]))
  } # end CU i
print(Sys.time() - start.time)

write.csv(tempThresh, file = "output/tempThresholds_FraserSEL.csv")

###############################################################################
# PLOTS
###############################################################################

# Plot map for each CU and lifestage
pdf(file = paste0("output/figures/Ts.percdays_FraserSEL.pdf"), width = 10, height = 5)

for(i in 1:length(cuids)){
  for(j in 1:length(stages)){
    
    # Subset temporal and spatial data
    cu_boundary.i <- cu_boundary[which(cu_boundary$FULL_CU_IN == cu_list$Full.CU.Index[which(cu_list$CUID == cuids[i])]),]
    
    zoi.i <- spawn_zoi[which(spawn_zoi$cuid == cuids[i]), ]
    
    mig_paths.i <- mig_paths[which(mig_paths$cuid == cuids[i]), ]
    
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
    
    # Highlight those cells
    grid_polys_incl <- st_as_sf(data.frame(
      V2 = rep(spat_grid$V2[incl], each = 5),
      rep(c("SW0", "NW", "NE", "SE", "SW1"), length(incl)),
      lon = c(rep(spat_grid$lon[incl], each = 5) + rep(c(- d/2, -d/2, d/2,  d/2, -d/2), length(incl))),
      lat = c(rep(spat_grid$lat[incl], each = 5) + rep(c(- d/2,  d/2, d/2, -d/2, -d/2), length(incl)))
    ), coords = c("lon", "lat"), crs = 4269) %>% 
      group_by(V2) %>%
      summarise(geometry = st_combine(geometry)) %>%
      st_cast("POLYGON") 
    
    
    # Set plot extent
    bounds <- st_bbox(grid_polys_incl)
    
    # Setup plot
    par(bg = NA, mfrow = c(2, 4), mar = c(0, 0, 0, 0), oma = c(2,4,4,2))
    
    for(k in 1:2){
      for(p in 1:4){
        
        plot(st_geometry(BC), 
             border = NA, 
             col = NA, 
             axes = FALSE, 
             las = 1, 
             xlim = bounds[c(1,3)], 
             ylim = bounds[c(2,4)], 
             bty = "o")
        u <- par('usr')
        segments(x0 = u[c(1, 1, 1, 2)], x1 = u[c(2, 2, 1, 2)], y0 = u[c(1, 2, 1, 1)], y1 = u[c(1, 2, 2, 2)], xpd = NA)
        
        if(k == 1) mtext(side = 3, line = 0, c("Historical", "Early-century", "Mid-century", "Late-century")[p])
        if(p == 1) mtext(side = 2, line = 2, c("Optimal", "Critical")[k])
        
        col_levels <- seq(0, 0.9, 0.1)
        col_palette <- colorRampPalette(colors = cols[c(5,3,1)])(n = length(col_levels) + 1)
        col_polys <- col_palette[findInterval(fw_exposure.spat[[i,j]][c("temp_opt", "temp_crit")[k], p, ], col_levels)]
        plot(grid_polys_incl, col = col_polys, border = col_polys, add = TRUE)
        
        # plot(grid_polys, add = TRUE, border = grey(0.8), col = NA, lwd = 0.5)
        
        # Plot coastline, lakes, and rivers (coarse)
        plot(BC, add = TRUE, col = NA, border = 1)
        plot(st_geometry(lakes0), border = 1, col = NA, add = TRUE, lwd = 0.7)
        plot(st_geometry(rivers0), col = 1, add = TRUE, lwd = 0.7)
        
        # legend("topright", fill = col_palette, legend = col_levels, bg = "white", title = "% days")
        
      }}
    mtext(side = 3, line = 1.5, paste(cus_to_keep$culabel[i], c("Adult freshwater migration", "Spawning", "Eggs/alevin", "Freshwater juvenile rearing")[j], sep = "-"), outer = TRUE)
  
}}
dev.off()

#-----------------------------------------------------------------------------
# Illustrate CU-level temperature thresholds
#-----------------------------------------------------------------------------

n.cu <- length(cuid)

#-----------------------------------------------------------------------------
# Comparing CUs
#-----------------------------------------------------------------------------
library(gplots)
# Historical

periodPch <- c(21, 25, 23, 24)
quartz(width = 8, height = 5, pointsize = 10)
p <- 1
par(mar = c(18, 4, 1, 1))
plotCI(1:length(cuids)-0.1, streamTemp_dd[, p, 1], li = streamTemp_dd[, p, 2], ui = streamTemp_dd[, p, 3], pch = periodPch[p], pt.bg = grey(0.6), col = grey(0.6), gap = 0, ylim = range(streamTemp_dd), xaxt = "n", xlim = c(0.5, 12.5), sfrac = 0, lwd = 1.5, las = 1, ylab = "", xlab = "", bty = "l")
axis(side = 1, at = 1:length(cuids), cus_to_keep$culabel, las = 2)

plotCI(1:length(cuids) + 0.5 - (0.05*p + 0.05), streamTemp_dd2[, p, 1], li = streamTemp_dd2[, p, 2], ui = streamTemp_dd2[, p, 3], pch = periodPch[p], pt.bg = popCol[cus_to_keep$population], col = popCol[cus_to_keep$population], gap = 0, sfrac = 0, lwd = 1.5, add = TRUE)

for(p in 2:4){
  plotCI(1:length(cuids) - (0.05*p + 0.05), streamTemp_dd[, p, 1], li = streamTemp_dd[, p, 2], ui = streamTemp_dd[, p, 3], pch = periodPch[p], pt.bg = grey(0.6), col = grey(0.6), gap = 0, sfrac = 0, lwd = 1.5, add = TRUE)
  plotCI(1:length(cuids) + 0.5 - (0.05*p + 0.05), streamTemp_dd2[, p, 1], li = streamTemp_dd2[, p, 2], ui = streamTemp_dd2[, p, 3], pch = periodPch[p], col = popCol[cus_to_keep$population], pt.bg = popCol[cus_to_keep$population],  gap = 0, sfrac = 0, lwd = 1.5, add = TRUE)
}
  
plot(1,1, "n")
legend("center", pch = periodPch, col = grey(0.6), pt.bg = grey(0.6), c("Historical", "Early", "Mid", "Late"))

Z_streamTemp_dd <- (streamTemp_dd - mean(streamTemp_dd))/sd(streamTemp_dd)
Z_streamTemp_dd2 <- (streamTemp_dd2 - mean(streamTemp_dd2[which(streamTemp_dd2>0)]))/sd(streamTemp_dd2[which(streamTemp_dd2>0)])

# streamTemp_dd <- Z_streamTemp_dd
# streamTemp_dd2 <- Z_streamTemp_dd2
# 
# streamTemp_dd0 <- streamTemp_dd
# streamTemp_dd20 <- streamTemp_dd2
# 
# streamTemp_dd <- streamTemp_dd0
# streamTemp_dd2<- streamTemp_dd20


#-----------------------------------------------------------------------------
# OLD
#-----------------------------------------------------------------------------

if(3 == 2){
  
  # MAP
  par(bg = NA)
  plot(st_geometry(zoi.i), border = cols[3], col = paste0(cols[3], 50), lwd = 1.5, axes= TRUE, las = 1)
  # plot(st_geometry(cu_boundary.i), border = cols[4], col = paste0(cols[4], 50), lwd = 1.5, add = TRUE)
  plot(st_geometry(cu_boundary.i), col = NA, lwd = 0.8, add = TRUE)
  
  u <- par('usr')
  bounds <- c(xmin = u[1], ymin = u[3], xmax = u[2], ymax = u[4])
  lakes.i <- lakes %>% st_crop(bounds)
  rivers.i <- rivers %>% st_crop(bounds)
  
  plot(st_geometry(lakes.i), border = cols[5], col = paste0(cols[5], 50), add = TRUE)
  plot(st_geometry(rivers.i), col = cols[5], add = TRUE)
  # plot(st_geometry(spawning_points), col = cols[1], pch = 19, add = TRUE, cex = 0.7)
  # plot(st_geometry(spawning_lines), col = cols[1], add = TRUE)
  
  # Add grid to map
  # plot(st_geometry(grid_polys), border = "#00000030", col = "#FFFFFF50", lwd = 0.5, add = TRUE)
  # plot(st_geometry(grid_polys_incl), border = "#00000080", add = TRUE)
  
  plot(st_geometry(grid_polys_incl)[[18]], col = "#00000080", add = TRUE)
  
  plot(st_geometry(st_multipoint(grid_points[[1]][-incl,])), add = TRUE, cex = 0.6)
  plot(st_geometry(st_multipoint(grid_points[[1]][incl,])), add = TRUE, cex = 0.6, pch = 19)
  
  #-----------------------------------------------------------------------------
  # Plot entire timeseries for one grid cell
  #-----------------------------------------------------------------------------
  
  par(mar = c(4,4,1,1))
  plot(date[seq(1, length(date), 7)], Ts[incl[18], seq(1, length(date), 7)], "l", xlab = "", ylab = "Stream temperature (˚C)", las = 1, xaxs = "i", yaxs = "i", bty = "l")
  abline(v = date[unlist(c(periods[1:4,4:5]))], col = cols[2], lwd = 2)
  
  plot(date[seq(1, length(date), 7)], Qs[incl[18], seq(1, length(date), 7)], "l", xlab = "", ylab = "Stream flow (cm/s)", las = 1, xaxs = "i", yaxs = "i", bty = "l")
  
  # Period - mid
  plot(date[seq(periods[3, 'start_num'], periods[3, 'end_num'], 1)], Ts[incl[18], seq(periods[3, 'start_num'], periods[3, 'end_num'], 1)], "l", xlab = "", ylab = "Stream temperature (˚C)", las = 1, xaxs = "i", yaxs = "i", bty = "l")
  plot(date[seq(periods[3, 'start_num'], periods[3, 'end_num'], 1)], Qs[incl[18], seq(periods[3, 'start_num'], periods[3, 'end_num'], 1)], "l", xlab = "", ylab = "Stream flow (cm/s)", las = 1, xaxs = "i", yaxs = "i", bty = "l")
  
  # Single year
  plot(date[which(years == 2060)], Ts[incl[18], which(years == 2060)], "l", xlab = "", ylab = "Stream temperature (˚C)", las = 1, xaxs = "i", yaxs = "i", bty = "l")
  lines(date[which(years == 2060)], weeklyMax(Ts[incl[18], which(years == 2060)]), lwd = 2, xpd = NA)
  # abline(h = topt[topt$Species == "Sockeye" & topt$Life.stage == "Spawning", c('min', 'max')], col = 2)
  
  plot(date[which(years == 1970)], Ts[incl[18], which(years == 1970)], "l", xlab = "", ylab = "Stream temperature (˚C)", las = 1, xaxs = "i", yaxs = "i", bty = "l")
  lines(date[which(years == 1970)], weeklyMax(Ts[incl[18], which(years == 1970)]), lwd = 2, xpd = NA)
  abline(h = 18, col = 2)
  
  
  abline(v = date[DOY %in% timing.i], col = cols[2])
  plot(date[which(years == 2060)], Qs[incl[18], which(years == 2060)], "l", xlab = "", ylab = "Stream flow (cm/s)", las = 1, xaxs = "i", yaxs = "i", bty = "l")
  
  # Plot with just spawn window
  plot(date[which(DOY %in% c(as.numeric(timing.i[1]):as.numeric(timing.i[2])))], Ts[incl[18], which(DOY %in% c(as.numeric(timing.i[1]):as.numeric(timing.i[2])))], "n", pch = 19, cex = 0.3, xlab = "", ylab = "Stream temperature (˚C)", las = 1, xaxs = "i", yaxs = "i", bty = "l")
  abline(v = date[unlist(c(periods[1:4,4:5]))], col = cols[2], lwd = 2)
  # abline(h = 18, col = 2)
  # for(j in 1:31){
  #   points(date[which(DOY %in% c(as.numeric(timing.i[1]):as.numeric(timing.i[2])))], Ts[incl[j], which(DOY %in% c(as.numeric(timing.i[1]):as.numeric(timing.i[2])))], pch = 19, cex = 0.3, col = "#00000030")
  # }
  for(j in 1:31){
    segments(x0 = date[periods$start_num], y0 = Ts.weeklyMax_avg[, j], x1 = date[periods$end_num], y1 = Ts.weeklyMax_avg[, j], col = cols[4], lwd = 1)
  }
  segments(x0 = date[periods$start_num], y0 = Ts.weeklyMax_avg[, 18], x1 = date[periods$end_num], y1 = Ts.weeklyMax_avg[, 18], col = cols[4], lwd = 3)
  
  # Plot avg through period
  plot(1:4, apply(Ts.dd, 1, mean), col = cols[4], pch = 19, ylim = range(Ts.dd), ylab = "Degree days above\n18˚C during spawning", las = 1, bty = "l", xaxt = "n", xlab = "", xlim = c(0.5,4.5), cex = 1.5)
  axis(side = 1, at = c(1:4), labels = c("Hist", "Early", "Mid", "Late"))
  segments(x0 = 1:4, x1 = 1:4, y0 = apply(Ts.dd, 1, min), y1 = apply(Ts.dd, 1, max), col = paste0(cols[4], 30), lwd = 5)
  for(j in 1:31){
    points(1:4, Ts.dd[,j], col = paste0(cols[4], 50), cex = 0.7)
  }
  
  # Change in DD
  delta_Ts.dd <- Ts.dd[2:4, ] - rbind(Ts.dd[1, ], Ts.dd[1, ], Ts.dd[1, ])
  plot(1:3, apply(delta_Ts.dd, 1, mean), col = cols[4], pch = 19, ylim = range(delta_Ts.dd), ylab = "Change in DD>18˚C \nduring spawning", las = 1, bty = "l", xaxt = "n", xlab = "", xlim = c(0.5,3.5), cex = 1.5)
  axis(side = 1, at = c(1:3), labels = c("Early", "Mid", "Late"))
  segments(x0 = 1:3, x1 = 1:3, y0 = apply(delta_Ts.dd, 1, min), y1 = apply(delta_Ts.dd, 1, max), col = paste0(cols[4], 30), lwd = 5)
  for(j in 1:31){
    points(1:3, delta_Ts.dd[,j], col = paste0(cols[4], 50), cex = 0.7)
  }
  
}
