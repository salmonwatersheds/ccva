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

source("freshwater/code/freshwater-functions.R")

# Conservation Units
cu_list <- read.csv("data/conservationunits_decoder.csv") %>% 
  subset(region == "Fraser" & (cu_type != "Extinct" | is.na(cu_type)))

# Check that CUs equal what's in incl.stages
cu_list$pooledcuid %in% as.numeric(dimnames(incl.stages)[[1]])
as.numeric(dimnames(incl.stages)[[1]]) %in% cu_list$pooledcuid
dim(cu_list)[1] == dim(incl.stages)[1]

# Create vector of cuid
cuid <- as.numeric(dimnames(incl.stages)[[1]])
n.CUs <- length(cuid)

# Order CU list to match incl.stages order (CK, CM, CO, PKE, PKO, SEL, SER, SH)
cu_list <- cu_list[match(cuid, cu_list$cuid), ]

# Create pooled species field (for matching to temp max)
cu_list$species_pooled <- cu_list$species_name
cu_list$species_pooled[cu_list$species_name %in% c("Lake sockeye", "River sockeye")] <- "Sockeye"
cu_list$species_pooled[cu_list$species_name %in% c("Pink (odd)", "Pink (even)")] <- "Pink"

# Spatial grid
grid_points <- read.csv("freshwater/data/processed-data/PCIC-grid-points_fraser.csv")

# Grid cells corresponding to each CU and stage for Fraser SEL
stages <- c("adult_migration", "spawning", "eggs_alevin", "fw_rearing")
incl.stages <- readRDS("freshwater/output/PCIC_incl.rds")
incl.stages_narm <- readRDS("freshwater/output/PCIC_incl_narm.rds")
dim(incl.stages)

# Subset Bowron-ES
i <- grep("Bowron", cu_list$cu_name_pse) #47
i <- grep("Upper Fraser", cu_list$cu_name_pse) #11
i <- which(grepl("Boundary Bay", cu_list$cu_name_pse) & cu_list$species_abbr == "CK") #1

incl <- incl.stages_narm[[i, which(stages == "adult_migration")]]# incl.stages[[i, which(stages == "adult_migration")]]


rcp <- 45
m <- 2

GCM_var <- readRDS(paste0("freshwater/data/processed-data/PCIC_", gcms$modelName[m], "_rcp", rcp,  "_processed.rds")) 

Ts.weeklyMax <- GCM_var[[1]]
Qs.weeklyMean <- GCM_var[[2]]
grid.ref <- as.numeric(dimnames(Ts.weeklyMax)[[1]])


start.ind <- which(period == c("hist", "early", "mid", "late")[p] & years == period.years[y, p] & DOY == timing.ij$start)

# Choose grid cell
mad <- mean(apply(Qs.weeklyMean[match(incl, grid.ref), which(period == "hist")], 2, median, na.rm = TRUE), na.rm = TRUE)


# Historical flow: median across grid cells and years
Q_hist0 <- apply(Qs.weeklyMean[match(incl, grid.ref), which(period == "hist")], 2, median, na.rm = TRUE)
DOY_hist <- strftime(as.Date(names(Qs.weeklyMean[1, which(period == "hist")])), format = "%b-%d")
Q_hist <- cbind(
  median = tapply(Q_hist0, DOY_hist, median),
  X25 = tapply(Q_hist0, DOY_hist, quantile, 0.25),
  X75 = tapply(Q_hist0, DOY_hist, quantile, 0.75))
DOY_plot <- as.Date(rownames(Q_hist), format = "%b-%d")
Q_hist <- Q_hist[order(DOY_plot),]

# Mid-century flow
Q_mid0 <- apply(Qs.weeklyMean[match(incl, grid.ref), which(period == "mid")], 2, median, na.rm = TRUE)
DOY_mid <- strftime(as.Date(names(Qs.weeklyMean[1, which(period == "mid")])), format = "%b-%d")
Q_mid <- cbind(
  median = tapply(Q_mid0, DOY_mid, median),
  X25 = tapply(Q_mid0, DOY_mid, quantile, 0.25),
  X75 = tapply(Q_mid0, DOY_mid, quantile, 0.75))
Q_mid <- Q_mid[order(DOY_plot),]


DOY_plot <- sort(DOY_plot)
DOY_tck <- as.Date(c(paste(2024, c(1, 4, 7, 10), "01", sep = "-"), "2024-12-31"))

quartz(width = 6, height = 4, pointsize = 10)
par(mar = c(4,5,3,1))
plot(DOY_plot, Q_hist[,'median'], "n", xaxt = "n", xlab = "", las = 1, ylab = expression(paste("Daily mean stream flow (", m^3, "/s)")), yaxs = "i", log = "y", xaxs = "i", ylim = c(0.1,30))

# add spawning period
timing <- read.csv("data/timing/timing-fraser.csv")
timing.ij <- timing[timing$cuid == cuid & timing$stage == "spawning", c('start', 'duration')]
x <- DOY_plot[as.numeric(c(timing.ij[1], timing.ij[1] + timing.ij[2]))]
polygon(x = c(x, rev(x)),
        y = c(10^-10, 10^-10, 30, 30),
        col = grey(0.9), border = NA)
text(mean(x), 0.15, "Spawn\ntiming", col = grey(0.4))

# add incubation?
timing.ij <- timing[timing$cuid == cuid & timing$stage == "eggs_alevin", c('start', 'duration')]
x <- DOY_plot[as.numeric(c(timing.ij[1], 366, 1, timing.ij[2] - (366 - timing.ij[1])))]
polygon(x = x[c(1,2,2,1)],
        y = c(10^-10, 10^-10, 30, 30),
        col = "#00000010", border = NA)
polygon(x = x[c(3,4,4,3)],
        y = c(10^-10, 10^-10, 30, 30),
        col = "#00000010", border = NA)
text(as.Date("2024-11-10"), 0.12, "Incubation", col = grey(0.4))
text(as.Date("2024-03-10"), 0.12, "Incubation", col = grey(0.4))

axis(side =1, at = DOY_tck, labels = c("Jan", "Apr", "Jul", "Oct", "Jan"))
abline(v = as.Date(paste(2024, c(2:12), "01", sep = "-")), lty = 3, col = grey(0.6))

# White polygons for flow
polygon(x = c(DOY_plot, rev(DOY_plot)),
        y = c(Q_hist[,'X25'], rev(Q_hist[,'X75'])),
        col = "#FFFFFF90",
        border = NA)
polygon(x = c(DOY_plot, rev(DOY_plot)),
        y = c(Q_mid[,'X25'], rev(Q_mid[,'X75'])),
        col = "#FFFFFF90",
        border = NA)

# Actual flow data with colour
lines(DOY_plot, Q_hist[,'median'], col = cols[5], lwd = 2)
polygon(x = c(DOY_plot, rev(DOY_plot)),
        y = c(Q_hist[,'X25'], rev(Q_hist[,'X75'])),
        col = paste0(cols[5], 30),
        border = NA)
abline(h = 0.2*mad, lty = 2)

lines(DOY_plot, Q_mid[,'median'], col = cols[4], lwd = 2)
polygon(x = c(DOY_plot, rev(DOY_plot)),
        y = c(Q_mid[,'X25'], rev(Q_mid[,'X75'])),
        col = paste0(cols[4], 30),
        border = NA)

legend("topleft", lwd = c(1,2,2), lty = c(2,1,1), col = c(1, cols[5], cols[4]), legend = c("20% MAD threshold", "Historical", "Mid century"), xpd = NA, bty = "n")

#------------------------------------------------------------------------------
# bcfishpass

# streams <- st_read(paste0(dat_root, "bcfishpass.gdb"), layer = "streams") %>%
  # st_transform(crs = 4269)

# Major rivers
# fwa <- st_read(paste0(XDrive_root, "/5_DATA/Mapping/study areas/shapefiles/hydrography_l_rivers_v2.shp")) %>%
  # st_transform(crs = 4269) # too coarse
wb <- st_read(paste0(XDrive_root, "/5_DATA/Mapping/FWA_Lakes_Poly_GT500HA_PSERegions.shp")) %>%
  filter(Region == "Fraser") %>%
  st_transform(crs = 4269) # too coarse
cu_boundary.i <- cu_boundary[which(cu_boundary$FULL_CU_IN == cu_list$cu_index[i]),]


plot(st_buffer(st_geometry(cu_boundary.i), 50*100), col = NA, border = NA)
plot(st_geometry(cu_boundary.i), border = cols[2], col = paste0(cols[2], 20), lwd = 2, add = TRUE)
plot(st_geometry(streams %>% filter(stream_order >= 5)), add = TRUE, col = cols[5], lwd = 2)
plot(st_geometry(streams %>% filter(stream_order >= 3)), add = TRUE, col = cols[5], lwd = 1.2)
plot(st_geometry(wb), col = cols[5], add = TRUE, border = NA)
plot(st_geometry(streams %>% filter(model_spawning_sk == 1)), add = TRUE)
plot(st_geometry(grid_polys_all), border = "#00000020", add = TRUE)
plot(st_geometry(grid_polys_all %>% filter(id %in% incl)), col = "#00000020", add = TRUE, border = NA)

plot(st_geometry(streams %>% filter(stream_order >= 5)), col = cols[5], lwd = 0.8)
plot(st_geometry(shoreline), add = TRUE)
plot(st_geometry(wb), col = cols[5], add = TRUE, border = cols[5])
plot(st_geometry(streams %>% filter(stream_order >= 9)), col = cols[5], lwd = 2, add = TRUE)
plot(st_geometry(cu_boundary.i), border = cols[2], col = paste0(cols[2], 20), lwd = 2, add = TRUE)

########################################################################################################
# Thompson summer steelhead CU - incubation
########################################################################################################

i <- grep("Thompson Summer", cu_list$cu_name_pse) #56

incl <- incl.stages_narm[[i, which(stages == "eggs_alevin")]]

# Temperature
# Historical flow: median across grid cells and years
T_hist0 <- apply(Ts.weeklyMax[match(incl, grid.ref), which(period == "hist")], 2, median, na.rm = TRUE)
DOY_hist <- strftime(as.Date(names(Ts.weeklyMax[1, which(period == "hist")])), format = "%b-%d")
T_hist <- cbind(
  median = tapply(T_hist0, DOY_hist, median),
  X25 = tapply(T_hist0, DOY_hist, quantile, 0.25),
  X75 = tapply(T_hist0, DOY_hist, quantile, 0.75))
DOY_plot <- as.Date(rownames(T_hist), format = "%b-%d")
T_hist <- T_hist[order(DOY_plot),]

# Mid-century flow
T_mid0 <- apply(Ts.weeklyMax[match(incl, grid.ref), which(period == "mid")], 2, median, na.rm = TRUE)
DOY_mid <- strftime(as.Date(names(Ts.weeklyMax[1, which(period == "mid")])), format = "%b-%d")
T_mid <- cbind(
  median = tapply(T_mid0, DOY_mid, median),
  X25 = tapply(T_mid0, DOY_mid, quantile, 0.25),
  X75 = tapply(T_mid0, DOY_mid, quantile, 0.75))
T_mid <- T_mid[order(DOY_plot),]

DOY_plot <- sort(DOY_plot)
DOY_tck <- as.Date(c(paste(2024, c(1, 4, 7, 10), "01", sep = "-"), "2024-12-31"))

quartz(width = 6, height = 4, pointsize = 10)
par(mar = c(4,5,3,1))
plot(DOY_plot, T_hist[,'median'], "n", xaxt = "n", xlab = "", las = 1, ylab = expression(paste("7-day max. of daily mean temperature (", degree, "C)")), yaxs = "i", xaxs = "i", ylim = c(0, 25))

# add spawning period
timing <- read.csv("data/timing/timing-fraser.csv")
timing.ij <- timing[timing$cuid == cu_list$cuid[i] & timing$stage == "eggs_alevin", c('start', 'duration')]
x <- DOY_plot[as.numeric(c(timing.ij[1], timing.ij[1] + timing.ij[2]))]
polygon(x = c(x, rev(x)),
        y = c(10^-10, 10^-10, 30, 30),
        col = grey(0.9), border = NA)
text(mean(x), 1, "Incubation", col = grey(0.4))

axis(side =1, at = DOY_tck, labels = c("Jan", "Apr", "Jul", "Oct", "Jan"))
abline(v = as.Date(paste(2024, c(2:12), "01", sep = "-")), lty = 3, col = grey(0.6))

# White polygons for flow
polygon(x = c(DOY_plot, rev(DOY_plot)),
        y = c(T_hist[,'X25'], rev(T_hist[,'X75'])),
        col = "#FFFFFF90",
        border = NA)
polygon(x = c(DOY_plot, rev(DOY_plot)),
        y = c(T_mid[,'X25'], rev(T_mid[,'X75'])),
        col = "#FFFFFF90",
        border = NA)

# Actual flow data with colour
lines(DOY_plot, T_hist[,'median'], col = cols[5], lwd = 2)
polygon(x = c(DOY_plot, rev(DOY_plot)),
        y = c(T_hist[,'X25'], rev(T_hist[,'X75'])),
        col = paste0(cols[5], 30),
        border = NA)
abline(h = tmax[i, 'eggs_alevin'], lty = 2)

lines(DOY_plot, T_mid[,'median'], col = cols[4], lwd = 2)
polygon(x = c(DOY_plot, rev(DOY_plot)),
        y = c(T_mid[,'X25'], rev(T_mid[,'X75'])),
        col = paste0(cols[4], 30),
        border = NA)

legend("topleft", lwd = c(1,2,2), lty = c(2,1,1), col = c(1, cols[5], cols[4]), legend = c("Temperature threshold\n(steelhead incubation)", "Historical", "Mid century"), xpd = NA, bty = "n")

#------------------------------------------------------------------------------
# Maps
cu_boundary.i <- cu_boundary[which(cu_boundary$CUID == cu_list$cuid[i]),]


plot(st_buffer(st_geometry(cu_boundary.i), 50*100), col = NA, border = NA)
plot(st_geometry(cu_boundary.i), border = cols[2], col = paste0(cols[2], 20), lwd = 2, add = TRUE)
plot(st_geometry(streams %>% filter(stream_order >= 6)), add = TRUE, col = cols[5], lwd = 2)
plot(st_geometry(streams %>% filter(stream_order >= 4)), add = TRUE, col = cols[5], lwd = 1.2)
plot(st_geometry(wb), col = cols[5], add = TRUE, border = NA)
plot(st_geometry(streams %>% filter(model_spawning_st == 1)), add = TRUE)
plot(st_geometry(grid_polys_all), border = "#00000020", add = TRUE, lwd = 0.8)
plot(st_geometry(grid_polys_all %>% filter(id %in% incl)), col = "#00000030", add = TRUE, border = NA)

plot(st_geometry(streams %>% filter(stream_order >= 5)), col = cols[5], lwd = 0.8)
plot(st_geometry(shoreline), add = TRUE)
plot(st_geometry(wb), col = cols[5], add = TRUE, border = cols[5])
plot(st_geometry(streams %>% filter(stream_order >= 9)), col = cols[5], lwd = 2, add = TRUE)
plot(st_geometry(cu_boundary.i), border = cols[2], col = paste0(cols[2], 20), lwd = 2, add = TRUE)
########################################################################################################
# South Thompson-Adams Upper Chinook CU - spawning
########################################################################################################

i <- grep("South Thompson-Adams River Upper", cu_list$cu_name_pse) #19

incl <- incl.stages[[i, which(stages == "spawning")]]

# Temperature
# Historical flow: median across grid cells and years
T_hist0 <- apply(Ts.weeklyMax[match(incl, grid.ref), which(period == "hist")], 2, median, na.rm = TRUE)
DOY_hist <- strftime(as.Date(names(Ts.weeklyMax[1, which(period == "hist")])), format = "%b-%d")
T_hist <- cbind(
  median = tapply(T_hist0, DOY_hist, median),
  X25 = tapply(T_hist0, DOY_hist, quantile, 0.25),
  X75 = tapply(T_hist0, DOY_hist, quantile, 0.75))
DOY_plot <- as.Date(rownames(T_hist), format = "%b-%d")
T_hist <- T_hist[order(DOY_plot),]

# Mid-century flow
T_mid0 <- apply(Ts.weeklyMax[match(incl, grid.ref), which(period == "mid")], 2, median, na.rm = TRUE)
DOY_mid <- strftime(as.Date(names(Ts.weeklyMax[1, which(period == "mid")])), format = "%b-%d")
T_mid <- cbind(
  median = tapply(T_mid0, DOY_mid, median),
  X25 = tapply(T_mid0, DOY_mid, quantile, 0.25),
  X75 = tapply(T_mid0, DOY_mid, quantile, 0.75))
DOY_plot <- as.Date(rownames(T_mid), format = "%b-%d")
T_mid <- T_mid[order(DOY_plot),]

# Late-century temp
T_late0 <- apply(Ts.weeklyMax[match(incl, grid.ref), which(period == "late")], 2, median, na.rm = TRUE)
DOY_late <- strftime(as.Date(names(Ts.weeklyMax[1, which(period == "late")])), format = "%b-%d")
T_late <- cbind(
  median = tapply(T_late0, DOY_late, median),
  X25 = tapply(T_late0, DOY_late, quantile, 0.25),
  X75 = tapply(T_late0, DOY_late, quantile, 0.75))
DOY_plot <- as.Date(rownames(T_late), format = "%b-%d")
T_late <- T_late[order(DOY_plot),]

DOY_plot <- sort(DOY_plot)
DOY_tck <- as.Date(c(paste(2024, c(1, 4, 7, 10), "01", sep = "-"), "2024-12-31"))

quartz(width = 6, height = 4, pointsize = 10)
par(mar = c(4,5,3,1))
plot(DOY_plot, T_hist[,'median'], "n", xaxt = "n", xlab = "", las = 1, ylab = expression(paste("7-day max. of daily mean temperature (", degree, "C)")), yaxs = "i", xaxs = "i", ylim = c(0, 35))

# add spawning period
timing <- read.csv("data/timing/timing-fraser.csv")
timing.ij <- timing[timing$cuid == cu_list$cuid[i] & timing$stage == "spawning", c('start', 'duration')]
x <- DOY_plot[as.numeric(c(timing.ij[1], timing.ij[1] + timing.ij[2]))]
polygon(x = c(x, rev(x)),
        y = c(10^-10, 10^-10, 30, 30),
        col = grey(0.9), border = NA)
text(mean(x), 1.5, "Spawning", col = grey(0.4))

axis(side =1, at = DOY_tck, labels = c("Jan", "Apr", "Jul", "Oct", "Jan"))
abline(v = as.Date(paste(2024, c(2:12), "01", sep = "-")), lty = 3, col = grey(0.6))

# White polygons for flow
polygon(x = c(DOY_plot, rev(DOY_plot)),
        y = c(T_hist[,'X25'], rev(T_hist[,'X75'])),
        col = "#FFFFFF90",
        border = NA)
polygon(x = c(DOY_plot, rev(DOY_plot)),
        y = c(T_mid[,'X25'], rev(T_mid[,'X75'])),
        col = "#FFFFFF90",
        border = NA)

# Actual flow data with colour
lines(DOY_plot, T_hist[,'median'], col = cols[5], lwd = 2)
polygon(x = c(DOY_plot, rev(DOY_plot)),
        y = c(T_hist[,'X25'], rev(T_hist[,'X75'])),
        col = paste0(cols[5], 30),
        border = NA)
# abline(h = tmax[i, 'spawning'], lty = 2)

lines(DOY_plot, T_mid[,'median'], col = cols[4], lwd = 2)
polygon(x = c(DOY_plot, rev(DOY_plot)),
        y = c(T_mid[,'X25'], rev(T_mid[,'X75'])),
        col = paste0(cols[4], 30),
        border = NA)

lines(DOY_plot, T_late[,'median'], col = cols[1], lwd = 2)
polygon(x = c(DOY_plot, rev(DOY_plot)),
        y = c(T_late[,'X25'], rev(T_late[,'X75'])),
        col = paste0(cols[1], 30),
        border = NA)

legend("topleft", lwd = c(1,2,2), lty = c(2,1,1), col = c(1, cols[5], cols[4]), legend = c("Temperature threshold\n(Chinook spawning)", "Historical", "Mid century"), xpd = NA, bty = "n")

#------------------------------------------------------------------------------
# Maps
cu_boundary.i <- cu_boundary[which(cu_boundary$CUID == cu_list$cuid[i]),]


plot(st_buffer(st_geometry(cu_boundary.i), 50*100), col = NA, border = NA, xlim = c(-120.5, -118.5))
plot(st_geometry(cu_boundary.i), border = cols[2], col = paste0(cols[2], 20), lwd = 2, add = TRUE)
plot(st_geometry(streams %>% filter(stream_order >= 6)), add = TRUE, col = cols[5], lwd = 2)
plot(st_geometry(streams %>% filter(stream_order >= 4)), add = TRUE, col = cols[5], lwd = 1.2)
plot(st_geometry(wb), col = cols[5], add = TRUE, border = NA)
plot(st_geometry(streams %>% filter(model_spawning_ch == 1)), add = TRUE)
plot(st_geometry(grid_polys_all), border = "#00000020", add = TRUE, lwd = 0.8)
plot(st_geometry(grid_polys_all %>% filter(id %in% incl)), col = "#00000030", add = TRUE, border = NA)

plot(st_geometry(streams %>% filter(stream_order >= 5)), col = cols[5], lwd = 0.8)
plot(st_geometry(shoreline), add = TRUE)
plot(st_geometry(wb), col = cols[5], add = TRUE, border = cols[5])
plot(st_geometry(streams %>% filter(stream_order >= 9)), col = cols[5], lwd = 2, add = TRUE)
plot(st_geometry(cu_boundary.i), border = cols[2], col = paste0(cols[2], 20), lwd = 2, add = TRUE)

###############################################################################
# Boundary Bay CK low flows during adult migration...
###############################################################################


###############################################################################
# Across stages, which are changing?

out.dat <- read.csv("output/ignore/exposure-raw.csv") %>%
  filter(rcp == "rcp85")

range(out.dat$median)
range(out.dat$max)
range(out.dat$min)

stage.names <- c("Incubation", "Freshwater rearing", "Early marine", "Marine rearing", "Adult migration", 'Spawning')
var.names <- cbind(c("stream temeprature", "stream flow"), c("SST", "salinity"))

quartz(width = 7.1, height = 9.5, pointsize = 10)
par(mfrow = c(3, 4), oma = c(1,10,1,0))
for(j in 1:length(stages)){
  for(v in 1:2){
    if (v == 1){
      par(mar = c(4, 1, 3, 0))
    } else {
      par(mar = c(4, 0, 3, 1))
    }
    
    var <- list(c("stream_temp", "stream_flow")[v], c("SST", "SSS")[v])[[c(1, 1, 2, 2, 1, 1)[j]]]
    plot(c(1,4), c(1,7), "n", yaxs = "i", xaxs = "i", xlim = c(0.8, 4.2), ylim = c(1, 7), bty = "n", xaxt = "n", xlab = "", ylab = "", yaxt = "n")
    axis(side = 1, at = c(1:4), c("H", "E", "M", "L"))
    abline(h = 1:6, lty = 3, col = grey(0.8))
    abline(v = 1:4, lty = 3, col = grey(0.8))
    abline(v = 3, lty = 3)
    mtext(side = 3, var.names[v, c(1,1,2,2,1,1)[j]], line = 0.2)
    if(v == 1) text(4.2, 7.9, xpd = NA, stage.names[j], cex = 1.6)
    if (j %in% c(1,3,5) & v == 1){
      mtext(side = 2, line = 8, "Exposure")
      for(s in 1:6){
        axis(side = 2, at = 7 - c(s, s - 1), label = c(0,1), las = 1, line = rep(c(-0.5, 1.5), 3)[s], col = species_cols[s])
      text(rep(c(0.6, 0.1), 3)[s], 7 - s + 0.5, names(species_cols)[s], col = species_cols[s], xpd = NA, adj = 1, cex = 1.3)
    }}
    for(s in 1:6){
      spql <- list("CK", "CM", "CO" ,"PKO", c("SEL", "SER"), "SH")[[s]]
      out.dat.js <- out.dat %>% filter(species %in% spql, stage == stages[j], variable == var)
      for(i in unique(out.dat.js$cuid)){
        out.dat.jsi <- out.dat.js %>% filter(cuid == i)
        y <- cbind(median = tapply(out.dat.jsi$median, out.dat.jsi$period, median),
                   min = tapply(out.dat.jsi$median, out.dat.jsi$period, min),
                   max = tapply(out.dat.jsi$median, out.dat.jsi$period, max))
        lines(1:4, 7 - s + y[,1], col = species_cols[s])
        polygon(x = c(1:4, 4:1),
                y = c(7 - s + y[,2], rev(7 - s + y[,3])),
                col = paste0(species_cols[s], 20), border = NA)
      } 
    }# end s
    
  }} # end j, v

mtext(side = 1, line = 0, outer = TRUE, "Period")
mtext(side = 3, outer = TRUE, line = -1, "RCP 8.5", adj = 0)
###############################################################################
###############################################################################
# Function to plot temp and flow
###############################################################################
###############################################################################
