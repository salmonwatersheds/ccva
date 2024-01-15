###############################################################################
# Code to summarize NOAA-provided GCM output for projected changes in sea-
# surface salinity and sea surface temperature
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

source("marine/code/load-NOAA.R")
###############################################################################
# Load NOAA spatial grid
###############################################################################

# NOAA grid points
grid <- read.csv("marine/data/raw-data/NOAA/NOAA-grid-points_wRegion_incl.csv")
grid_coast400 <- read.csv("marine/data/raw-data/NOAA/NOAA-grid-points_wRegion_incl_coast.csv")
grid_coast1000 <- read.csv("marine/data/raw-data/NOAA/NOAA-grid-points_wRegion_incl_coast1000.csv")
#%>%
  # subset(hasData == 1)
# grid_points <- st_as_sf(NOAA_grid, coords = c("lon", "lat"), crs = 4269)

# What are the unique NOAA grid ids that we'd want to include
unique_id <- unique(c(grid$id, grid_coast400$id, grid_coast1000$id))


###############################################################################
# Calculate Z scores for early, mid, and late century for SST
###############################################################################

sst <- loadNOAA(variable = "SST", model = "Can-ESM2", spat_id = unique_id)
sst_ids <- as.numeric(dimnames(sst)[[1]])

dates <- as.Date(dimnames(sst)[[2]])
months <- as.numeric(strftime(dates, format = "%m"))
period_ind <- list(
  hist = which(dates >= as.Date("1970-01-01") & dates < as.Date("2000-01-01")),
  early = which(dates >= as.Date("2010-01-01") & dates < as.Date("2040-01-01")),
  mid = which(dates >= as.Date("2040-01-01") & dates < as.Date("2070-01-01")),
  late = which(dates >= as.Date("2070-01-01") & dates < as.Date("2100-01-01"))
)

# For each unique grid cell, calculate mean and sd over historical period
hist.mean <- apply(sst[, period_ind$hist], 1, mean, na.rm = TRUE)
hist.sd <- apply(sst[, period_ind$hist], 1, sd, na.rm = TRUE)

z <- array(NA, dim = c(dim(sst)[1], 12, 3))

for(j in 1:3){
  for(m in 1:12){
    ind <- period_ind[[j+1]][which(months[period_ind[[j+1]]] == m)]
    z[, m, j] <- (apply(sst[, ind], 1, mean) - hist.mean) / hist.sd
  }
}

###############################################################################
# Summarize metrics over different areas
###############################################################################

# List of ids for four scenarios:
# 1) Full 
# 2) Pacific proper
# 3) Bering Sea
# 4) Coastal 400
# 5) Coastal 1000

area_ids <- list(
  grid$id,
  grid$id[grid$region %in% c("Gulf of Alaska", "North Pacific Ocean")],
  grid$id[grid$region %in% c("Bering Sea")],
  grid_coast400$id,
  grid_coast1000$id
)

sum(area_ids[[1]] %in% sst_ids == FALSE) # There are some that must be DD

median_z <- array(NA, dim = c(5, 12, 3))
for(i in 1:5){
  for(j in 1:3){
  median_z[i, , j] <- apply(z[which(sst_ids %in% area_ids[[i]]), , j], 2, median)
  }}

quartz(width = 16, height = 7, pointsize = 16)
par(mfrow = c(2, 6), mar = c(1,1,3,1), oma= c(3,3,2,1))
for(m in 1:12){
  barplot(median_z[, m, 2], col = cols, main = month.abb[m], ylim = c(-0.6, 3))
  abline(h = 0)
  abline(h = c(-1.28, 1.28), lty = 2, lwd = 1.5)
}
plot(1,1,"n")
legend("center", fill = cols, legend = c("Full", "Pacific", "Bering", "Coastal 400", "Coastal 1000"))

# proportion of cells
prop_z <- array(NA, dim = c(5, 12, 3))
for(i in 1:5){
  for(j in 1:3){
    prop_z[i, , j] <- apply(z[which(sst_ids %in% area_ids[[i]]), , j], 2, function(x){length(which(x >= 1.28 | x <= -1.28))/length(x)})
  }}

quartz(width = 16, height = 7, pointsize = 16)
par(mfrow = c(2, 6), mar = c(1,1,3,1), oma= c(3,3,2,1))
for(m in 1:12){
  barplot(prop_z[, m, 2], col = cols, main = month.abb[m])
  abline(h = 0)
}
plot(1,1,"n")
legend("center", fill = cols, legend = c("Full", "Pacific", "Bering", "Coastal 400", "Coastal 1000"))

# Proportion of cells outside core range
coreRange <- data.frame(
  species = c("Chinook", "Chum", "Coho", "Pink", "Sockeye", "Steelhead"),
  lower = c(6.25, 3.25, 7, 3.8, 2.9, 6.5),
  upper = c(10.3, 10, 10.3, 8.75, 7.9, 10)
)

prop_core <- array(NA, dim = c(6, 5, 12, 3)) # species, region, month, period

  for(j in 1:3){
    for(m in 1:12){
      ind <- period_ind[[j+1]][which(months[period_ind[[j+1]]] == m)]
      mean.sst <- apply(sst[, ind], 1, mean)
      for(i in 1:5){
        mean.sst.i <- mean.sst[which(sst_ids %in% area_ids[[i]])] 
        for(s in 1:6){
          prop_core[s, i, m, j] <- length(which(mean.sst.i < coreRange$lower[s] | mean.sst.i > coreRange$upper[s]))/length(mean.sst.i)
        }}}}

quartz(width = 16, height = 7, pointsize = 16)
par(mfrow = c(2, 6), mar = c(1,1,3,1), oma= c(3,3,2,1))
for(s in 1:6){
  for(m in 1:12){
    barplot(prop_core[s, , m, 2], col = cols, main = month.abb[m])
    abline(h = 0)
  }
  mtext(side = 3, outer = TRUE, paste0(c("Chinook", "Chum", "Coho", "Pink", "Sockeye", "Steelhead")[s], " (", coreRange$lower[s], "-", coreRange$upper[s], "˚C)"))
}

# Average for stream-type Chinook
quartz(width = 4, height = 4, pointsize = 12)
barplot(apply(prop_core[1, , c(rep(1:12, 4), 8), 2], 1, mean), col = cols, las = 1, ylab = "Prop. outside core temp range")
barplot(apply(prop_core[1, , c(rep(1:12, 4), 8), 2], 1, sum), col = cols, las = 1, ylab = "Cum. outside core temp range", ylim = c(0, 40))
abline(h = 0)
