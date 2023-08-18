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
            dim = c(4, length(cuids), 3, 3),
            dimnames = list(
              c("temp_opt", "temp_crit", "lowflow_opt", "lowflow_crit"),
              cus_to_keep$culabel, 
              stage
              c("early", "mid", "late"), 
              c("median", "lower", "upper")))

#------------------------------------------------------------------------------
# For each CU
#------------------------------------------------------------------------------

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
      timing.ij <- timing[which(timing$cuid == cuids[i] & timing$stage == stages[j])[1], c("start", "end")] %>% 
      as.character() %>%
      as.Date(format = "%Y-%m-%d") %>% 
      strftime(format = "%j") %>% 
      as.numeric() %>%
      round()
  
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
      
  #-----------------------------------------------------------------------------
  # Calculate MAD over historical period
  #-----------------------------------------------------------------------------
  mad <- apply(Qs[incl, c(periods[periods$period == 'hist', 'start_num']:periods[periods$period == "hist", "end_num"])], 1, mean)
  
  #-----------------------------------------------------------------------------
  # Calculate rolling weekly max (Ts) or mean (Qs) over spawning period
  #-----------------------------------------------------------------------------
  
  # Select dates in spawn timing, extending earlier 7 days for calculating 7-day rolling avg
  periods.ind <- rbind(
    hist = which(DOY %in% c(c(timing.i$spawn_start - 7):timing.i$spawn_end) & years %in% c(1970:1999)),
    early = which(DOY %in% c(c(timing.i$spawn_start - 7):timing.i$spawn_end) & years %in% c(2010:2039)),
    mid = which(DOY %in% c(c(timing.i$spawn_start - 7):timing.i$spawn_end) & years %in% c(2040:2069)),
    late = which(DOY %in% c(c(timing.i$spawn_start - 7):timing.i$spawn_end) & years %in% c(2070:2099)))
  
  # Create arrays
  Ts.weeklyMax <- array(
    NA, 
    dim = c(4, length(incl), ncol(periods.ind) - 7), 
    dimnames = list(c("hist", "early", "mid", "late")))
  Qs.weeklyMean <- Ts.weeklyMax
  
  for(p in 1:4){ # For each period
    Ts.weeklyMax[p, , ] <- t(apply(Ts[incl, periods.ind[p,]], 1, weeklyStat, stat = "max"))[, 8:ncol(periods.ind)]
    Qs.weeklyMean[p, , ] <- t(apply(Qs[incl, periods.ind[p,]], 1, weeklyStat, stat = "mean"))[, 8:ncol(periods.ind)]
  }
  
  # Checks in case -Inf (not getting this any more?)
  Ts.weeklyMax[which(Ts.weeklyMax == -Inf)] <- NA
  Qs.weeklyMean[which(Qs.weeklyMean == -Inf)] <- NA
  
  # Average weekly max for each period
  # Calculate degree days (dd) over thresh
  Ts.weeklyMax_avg <- matrix(NA, nrow = 4, ncol = length(incl))
  Ts.dd <- matrix(NA, nrow = 4, ncol = length(incl))
  Ts.dd2 <- matrix(NA, nrow = 4, ncol = length(incl))
  for(p in 1:4){
    Ts.weeklyMax_avg[p, ] <- apply(Ts.weeklyMax[p, , ], 1, mean, na.rm = TRUE)
    Ts.dd[p, ] <- apply(Ts.weeklyMax[p, , ], 1, degdays, thresh = topt)
    Ts.dd2[p, ] <- apply(Ts.weeklyMax[p, , ], 1, degdays, thresh = as.numeric(topt2[topt2$population == cus_to_keep$population[i], c("min", "max")]))
  }
  
  # Summarize for CU
  streamTemp_dd[i, 1:4, 1] <- apply(Ts.dd, 1, median)
  streamTemp_dd[i, 1:4, 2] <- apply(Ts.dd, 1, min)
  streamTemp_dd[i, 1:4, 3] <- apply(Ts.dd, 1, max)
  
  streamTemp_dd2[i, 1:4, 1] <- apply(Ts.dd2, 1, median)
  streamTemp_dd2[i, 1:4, 2] <- apply(Ts.dd2, 1, min)
  streamTemp_dd2[i, 1:4, 3] <- apply(Ts.dd2, 1, max)
  
  #-----------------------------------------------------------------------------
  # Flow
  #-----------------------------------------------------------------------------
  
  # Calculate degree days (dd) over thresh
  Qs.cmsd <- matrix(NA, nrow = 4, ncol = length(incl))
  for(p in 1:4){
    for(j in 1:length(incl)){
      Qs.cmsd[p, j] <- degdays(Qs.weeklyMean[p, j, ], thresh = c(0, 0.2*mad[j]))
    }
  }
  
  
  # Summarize for CU
  lowFlow_cmsd[i, 1:4, 1] <- apply(Qs.cmsd, 1, median)
  lowFlow_cmsd[i, 1:4, 2] <- apply(Qs.cmsd, 1, min)
  lowFlow_cmsd[i, 1:4, 3] <- apply(Qs.cmsd, 1, max)
  
  }


###############################################################################
# PLOTS
###############################################################################

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