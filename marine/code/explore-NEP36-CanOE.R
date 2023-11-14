###############################################################################
# Code to explore downscaled marine climate data for coastal BC from
# Holdsworth et al. (2021), available at 
# https://open.canada.ca/data/en/dataset/5551969d-f94c-406b-b849-50b49d32290f 
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
library(PNWColors)

###############################################################################
# Import NetCDF file and explore
###############################################################################

# Open netcdf file
z <- nc_open("marine/data/raw-data/NEP36-CanOE/NEP36-CanOE_temp_RCP45_2046-2065_monthly.nc")

# Explore dataset. What are the variables?
print(paste("File",z$filename,"contains",z$nvars,"variables"))
for( i in 1:z$nvars ) {
  v <- z$var[[i]]
  print(paste("Here is information on variable number",i))
  print(paste("   Name: ",v$name))
  print(paste("   Units:",v$units))
  print(paste("   Missing value:",v$missval))
  print(paste("   # dimensions :",v$ndims))
  print(paste("   Variable size:",v$varsize))
}

# Extract metadata
sink(paste0("marine/data/raw-data/NEP36-CanOE/NEP36-CanOE_temp_RCP45_2046-2065_monthly", "_metadata.txt"))
print(z)
sink()


# Has temperature at depth; just use 1 m depth
depth <- ncvar_get(z, "deptht_bounds")
# First depth element is surface 0 to 1 m: use this.

temp.all <- ncvar_get(z, "temp")
dim(temp.all) # dimensions [x,y,deptht,t]
lat <- ncvar_get(z, "nav_lat")
lon <- ncvar_get(z, "nav_lon")
dim(lon)


sst <- ncvar_get(z, "temp")[, , 1, ]

# Map January sst
tempCols <- pnw_palette("Bay", n = 30) 
tempLevels <- seq(min(sst, na.rm = TRUE), max(sst, na.rm = TRUE), length.out = 29)
plot(lon[1,], lat[1,], pch = 19, col = tempCols[findInterval(sst[1, , 1], tempLevels)])

# Create grid points spatial
sst_points <- data.frame(
  id = c(1:(dim(sst)[1] * dim(sst)[2])),
  x = c(matrix(rep(1:715, 1021))),
  y = c(matrix(rep(1:1021, each = 715))),
  lon = c(lon),
  lat = c(lat)) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4269)
  
plot(st_geometry(sst_points), col = tempCols[findInterval(sst[1, , 1], tempLevels)], cex = 0.5)

#------------------------------------------------------------------------------
# Import future and historical sst
#------------------------------------------------------------------------------

# Open netcdf file
z.future <- nc_open("marine/data/raw-data/NEP36-CanOE/NEP36-CanOE_temp_RCP45_2046-2065_monthly.nc")
z.hist <- nc_open("marine/data/raw-data/NEP36-CanOE/NEP36-CanOE_temp_historical_1986-2005_monthly.nc")

sst.future <- ncvar_get(z.future, "temp")[, , 1, ]
sst.hist <- ncvar_get(z.hist, "temp")[, , 1, ]

#------------------------------------------------------------------------------
# Read in early marine distribution
#------------------------------------------------------------------------------
shoreline <- st_read(dsn = "freshwater/data/spatial/layers/GSHHS_i_L1.shp")

region <- "Fraser" # Choose between Fraser, Skeena, Nass
# Load early marine spatial zone (em_spat)
em_spat <- st_read(dsn = paste0("marine/data/spatial/early-marine-shapefiles/", region, "/", region, "_early_marine.shp")) %>% st_transform(crs = 4269)

plot(st_geometry(em_spat), col = grey(0.8), border = NA)
plot(st_geometry(shoreline), add = TRUE)

plot(st_geometry(em_spat), col = grey(0.8), border = NA, xlim = c(-123.8, -123), ylim = c(48, 48.2))
plot(st_geometry(shoreline), add = TRUE)
plot(st_geometry(sst_points), pch = 19, cex = 0.3, col = cols[1], add = TRUE)     

# Subset sst points in early marine spatial zone (em_spat)
intrscts <- st_intersects(sst_points, em_spat, sparse = FALSE)
incl.em <- which(c(intrscts) == TRUE)

head(sst_points[incl.em,])

# Calculate historical and future
sst <- array(NA, dim = c(length(incl.em), 12, 2))
for(j in 1:12){
  sst[, j, 1] <- sst.hist[cbind(sst_points$x[incl.em], sst_points$y[incl.em], j)]
  sst[, j, 2] <- sst.future[cbind(sst_points$x[incl.em], sst_points$y[incl.em], j)]
}
  
#------------------------------------------------------------------------------
# Read in early marine timing
#------------------------------------------------------------------------------
timing <- read.csv("data/timing/3Life_cycle_timing_by_CU.csv")

# Extract ocean entry
timing$oe <- apply(timing[, c('oe_start', "oe_end")], 1, mean, na.rm = TRUE)
oe <- timing$oe_end[which(timing$region == "fraser")]

# Ocean entry month
oe_month <- as.numeric(strftime(as.Date(paste(2001, oe, sep = "-"), format = "%Y-%j"), format = "%m"))
# Some variation in month of ocean entry, but most are May. Work with that for now

barplot(tapply(oe_month, oe_month, length), names = c("Apr", "May", "Jun", "Jul"), ylab = "Number of CUs", las = 1, xlab = "Ocean entry month")
abline(h = 0)

data.frame(
  region = unique(timing$region),
  feb_start = NA,
  mar_start = NA,
  apr_start = NA,
  may_start = NA,
  jun_start = NA,
  apr_end = NA,
  may_end,
  jun_end,
  jul_end,
  aug_end,
)

xDate <- as.Date(paste("2001", rep(c(2:6), each = 2), rep(c(1, 15), 5), sep = "-"))
xDOY <- as.numeric(strftime(xDate, format = "%j"))
xLabels <- strftime(xDate, format = "%d-%b") 
xBreaks <- as.numeric(strftime(as.Date(paste("2001", rep(c(2:6), each = 4), rep(c(1, 8, 15, 22), 5), sep = "-")), format = "%j"))

quartz(width = 8, height = 3.6, pointsize = 10)
hist(timing$oe_start, freq = TRUE, breaks = xBreaks, main = "Start of ocean entry", xaxt = "n", xlab = NA, col = cols[as.numeric(strftime(as.Date(paste("2001", rep(c(2:6), each = 4), rep(c(1, 8, 15, 22), 5), sep = "-")), format = "%m")) - 1], las = 1, yaxs = "i", ylab = "Number of Conservation Units")
axis(side = 1, at = xDOY, labels = xLabels)
legend("topleft", fill = cols[1:5], legend = c("Feb" ,"Mar", "Apr", "May", "Jun"), bty = "n")


# Break down by species
SQ <- unique(timing$species)
timing$species_pooled <- timing$species
timing$species_pooled[timing$species %in% c("CK")] <- "Chinook"
timing$species_pooled[timing$species %in% c("CM")] <- "Chum"
timing$species_pooled[timing$species %in% c("CO")] <- "Coho"
timing$species_pooled[timing$species %in% c("PKE", "PKO")] <- "Pink"
timing$species_pooled[timing$species %in% c("SEL", "SER")] <- "Sockeye"
timing$species_pooled[timing$species %in% c("SH")] <- "Steelhead"

quartz(width = 8, height = 8, pointsize = 10)
par(mfrow = c(6,1), mar = c(1, 4, 1, 1), oma = c(3, 0, 0, 0))
for(s in 1:6){
  hist(timing$oe_start[timing$species_pooled == unique(timing$species_pooled)[s]], freq = TRUE, breaks = xBreaks, main = "", xaxt = "n", xlab = NA, col = cols[as.numeric(strftime(as.Date(paste("2001", rep(c(2:6), each = 4), rep(c(1, 8, 15, 22), 5), sep = "-")), format = "%m")) - 1], las = 1, yaxs = "i", ylab = "Number of Conservation Units")
  mtext(side =3, line = -1, adj = 0, paste0("   ", unique(timing$species_pooled)[s]))
  # legend("topleft", fill = cols[1:5], legend = c("Feb" ,"Mar", "Apr", "May", "Jun"), bty = "n")
  axis(side = 1, at = xDOY, labels = FALSE)
  
}
axis(side = 1, at = xDOY, labels = xLabels)

#------------------------------------------------------------------------------
# Plot month by month Diff
#------------------------------------------------------------------------------

# Calculate z score of (future - mean(hist))/sd(hist)
z_score <- (sst[, , 2] - apply(sst[, , 1], 1, mean))/apply(sst[, , 1], 1, sd) 

quartz(width = 6, height = 4)
par(mar = c(4,4,2,1))
hist(sst[, 5, 1], breaks = seq(0, 30, 0.2), col = paste0(cols[5], 50), border = NA, xlim = c(6, 20), xlab = "Sea surface temperature (˚C)", main = "", yaxs = "i")
hist(sst[, 5, 2], breaks = seq(0, 30, 0.2), col = paste0(cols[1], 50), border = NA, add = TRUE)
legend("topright", fill = paste0(cols[c(5,1)], 50), legend = c("Historical (1986-2005)", "Future (2046-2065"), border = NA, bty = "n")


hist(z_score[, 7], breaks = seq(-3, 7, 0.1), xlim = c(-2, 3), col = c(rep(cols[5], 30), rep(cols[1], 71)), xlab = "z-score", main = "July", yaxs = "i")

zCols <- pnw_palette("Bay", n = 30) 
zLevels <- seq(-2, 3, length.out = 29)


plot(st_geometry(em_spat), col = NA, border = NA)
plot(st_geometry(sst_points[incl.em, ]), col = zCols[findInterval(z_score[, 10], zLevels)], pch = 19, cex = 0.1, add = TRUE)
plot(st_geometry(shoreline), add = TRUE)

plot(rep(1, 30), 1:30, pch = 19, col = zCols, cex = 2)
text(0.98, seq(1.5, 29.5, 1), round(zLevels, 2), pos = 2)
segments(x0 = 0.98, x1 = 1, y0 = seq(0.5, 29.5, 1))


z_out <- data.frame(
  oe_month = c("Apr", "May", "Jun", "Jul"),
  mean_z = c(mean(z_score[, 2:7], na.rm = TRUE), mean(z_score[, 3:8], na.rm = TRUE), mean(z_score[, 4:9], na.rm = TRUE), mean(z_score[, 5:10], na.rm = TRUE)),
  min_z = c(min(z_score[, 2:7], na.rm = TRUE), min(z_score[, 3:8], na.rm = TRUE), min(z_score[, 4:9], na.rm = TRUE), min(z_score[, 5:10], na.rm = TRUE)),
  max_z = c(max(z_score[, 2:7], na.rm = TRUE), max(z_score[, 3:8], na.rm = TRUE), max(z_score[, 4:9], na.rm = TRUE), max(z_score[, 5:10], na.rm = TRUE))
)

plot(1:4, z_out$mean_z, ylim = range(z_out[, 3:4]))
segments(x0 = 1:4, x1 = 1:4, y0 = z_out$min_z, y1 = z_out$max_z)

barplot(z_out[, 2], ylab = "mean z-score", xlab = "Ocean entry month", names.arg = z_out$oe_month, las = 1)


z_dens <- cbind(
  x = density(z_score[, 2:7], na.rm = TRUE, bw = 0.5, from = -3, to = 5, n = 100)$x,
  Apr = density(z_score[, 2:7], na.rm = TRUE, bw = 0.5, from = -3, to = 5, n = 100)$y,
  May = density(z_score[, 3:8], na.rm = TRUE, bw = 0.5, from = -3, to = 5, n = 100)$y,
  Jun = density(z_score[, 4:9], na.rm = TRUE, bw = 0.5, from = -3, to = 5, n = 100)$y,
  Jul = density(z_score[, 5:10], na.rm = TRUE, bw = 0.5, from = -3, to = 5, n = 100)$y
)

head(z_dens)

par(mar = c(3, 4, 4, 1))
plot(z_dens[, 'x'], z_dens[, "Apr"], "n",  xlab = "Z-score", ylab = "Density", bty = "l", las = 1)
abline(v = c(0.5, 2), lty = 2)
u <- par('usr')
for(i in 1:4){
  med <- median(z_score[, (c(2:5)[i]):(c(7:10)[i])], na.rm = TRUE)
  lines(z_dens[, 'x'], z_dens[, i+1], "l", col = cols[c(2,3,4,1)[i]], lwd = 1.5, xpd = NA)
  points(med, u[4] + (i - 1) *0.05*(u[4] - u[3]), col = cols[c(2,3,4,1)[i]], pch = 19, cex = 1.5, xpd = NA)
  segments(x0 = quantile(z_score[, (c(2:5)[i]):(c(7:10)[i])], 0.25, na.rm = TRUE),
           x1 = quantile(z_score[, (c(2:5)[i]):(c(7:10)[i])], 0.75, na.rm = TRUE),
           y0 = u[4] + (i - 1) *0.05*(u[4] - u[3]),
           y1 = u[4] + (i - 1) *0.05*(u[4] - u[3]),
           col = cols[c(2,3,4,1)[i]],
           xpd = NA)
}
legend("topleft", lwd = 1.5, col = cols[c(2,3,4,1)], c("Apr", "May", "Jun", "Jul"), bty = "n")


#------------------------------------------------------------------------------
# Compare to rough NOAA
#------------------------------------------------------------------------------

fraser_temps <- readRDS("marine/data/raw-data/NOAA/fraser-summary.rds")
fraser_incl <- read.csv("marine/data/raw-data/NOAA/Fraser_early_marine_ids.csv")

noaa_points <- data.frame(
  id = fraser_incl$id,
  lon = fraser_incl$lon - 360,
  lat = fraser_incl$lat) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4269)

plot(st_geometry(noaa_points), col = zCols[findInterval(z_score2[, 5], zLevels)], pch = 19, cex = 2)
plot(st_geometry(shoreline), add = TRUE, lwd = 0.5)

# Calculate z score of (future - mean(hist))/sd(hist)
z_score2 <- (fraser_temps[, , 2] - apply(fraser_temps[, , 1], 1, mean))/apply(fraser_temps[, , 1], 1, sd) 

z_out2 <- data.frame(
  oe_month = c("Apr", "May", "Jun", "Jul"),
  mean_z = c(mean(z_score2[, 2:7], na.rm = TRUE), mean(z_score2[, 3:8], na.rm = TRUE), mean(z_score2[, 4:9], na.rm = TRUE), mean(z_score2[, 5:10], na.rm = TRUE)),
  min_z = c(min(z_score2[, 2:7], na.rm = TRUE), min(z_score2[, 3:8], na.rm = TRUE), min(z_score2[, 4:9], na.rm = TRUE), min(z_score2[, 5:10], na.rm = TRUE)),
  max_z = c(max(z_score2[, 2:7], na.rm = TRUE), max(z_score2[, 3:8], na.rm = TRUE), max(z_score2[, 4:9], na.rm = TRUE), max(z_score2[, 5:10], na.rm = TRUE))
)
barplot(rbind(z_out[,2], z_out2[, 2]), beside = TRUE, ylab = "mean z-score", xlab = "Ocean entry month", names.arg = z_out$oe_month, las = 1, col = cols[c(2,3)])
abline(h = 0)
legend("topleft", fill = cols[c(2,3)], legend = c("NEP36-CanOE downscaled", "CanESM2 raw"), bty = "n")

###############################################################################
# Summary of z-scores of temp and sal changes for Nass, Skeena, Fraser by OE months
###############################################################################

grid <- read.csv("marine/data/raw-data/NEP36-CanOE/NEP36-CanOE_grid.csv")
# 715 x 1021 grid points

grid_points <- st_as_sf(grid, coords = c("lon", "lat"), crs = 4269)

z_scores <- data.frame(
  variable = rep(c('temp', 'salt'), each = 3*5),
  region = rep(rep(c("Fraser", "Skeena", "Nass"), each = 5), 2),
  oe_start = rep(c(2:6), 3*2),
  z025 = NA,
  z_25 = NA,
  z_median = NA,
  z_75 = NA,
  z_975 = NA
)


#--------
# For both variables: Sea surface temperature (temp) and salinity (salt)
#----------

for(i in 1:2){
  
  # Open netcdf file
  root <- "marine/data/raw-data/NEP36-CanOE/"
  z.future <- nc_open(paste0(root, "NEP36-CanOE_", c("temp", "salt")[i], "_RCP45_2046-2065_monthly.nc"))
  z.hist <- nc_open(paste0(root, "NEP36-CanOE_", c("temp", "salt")[i], "_historical_1986-2005_monthly.nc"))
  
  depth <- 1 # Select surface salinity and temp
  
  # Extract relevant variables
  var.future <- ncvar_get(z.future, c("temp", "salt")[i])[, , depth, ]
  var.hist <- ncvar_get(z.hist, c("temp", "salt")[i])[, , depth, ]
  
  #----------
  # For each region
  #----------
  
  for(r in 1:3){
    region <- unique(z_score.sst$region)[r]
    
    # Find downscaled grid points corresponding to region
    # Load early marine shapefile
    em_spat <- st_read(dsn = paste0("marine/data/spatial/early-marine-shapefiles/", region, "/", region, "_early_marine.shp")) %>% st_transform(crs = 4269)
    
    # Subset sst points in early marine spatial zone (em_spat)
    intrscts <- st_intersects(grid_points, em_spat, sparse = FALSE)
    incl.em <- which(c(intrscts) == TRUE)
    
    # # Plot to check
    # plot(st_geometry(em_spat), col = grey(0.8), border = NA)
    # plot(st_geometry(shoreline), add = TRUE)
    # # plot(st_geometry(sst_points), pch = 19, cex = 0.3, col = cols[1], add = TRUE)  
    # plot(st_geometry(sst_points[incl.em, ]), pch = 19, cex = 0.05, col = paste0(cols[1], 10), add = TRUE) 
    
    # Calculate historical and future
    var.region <- array(NA, dim = c(length(incl.em), 12, 2))
    for(j in 1:12){
      var.region[, j, 1] <- var.hist[cbind(grid$x[incl.em], grid$y[incl.em], j)]
      var.region[, j, 2] <- var.future[cbind(grid$x[incl.em], grid$y[incl.em], j)]
    }
    
    # Calculate z score of (future - mean(hist))/sd(hist)
    z_score <- (var.region[, , 2] - apply(var.region[, , 1], 1, mean))/apply(var.region[, , 1], 1, sd) 
    
    # Store spatially explicit output
    write.csv(cbind(grid$id[incl.em], z_score), file = paste0("marine/output/NEP36-CanOE/fullZscores_NEP36-CanOE", "_", region, "_", c("temp", "salt")[i], ".csv"))
  
    # Summarize
    for(j in 1:5){ # For each start OE month from Feb - Jun
      z_scores[which(z_scores$variable ==  c("temp", "salt")[i] & z_scores$region == region & z_scores$oe_start == j + 1), c("z025", "z_25", "z_75", "z_975")] <- quantile(z_score[, cbind(c(11:12,1:4), c(12,1:5), c(1:6), c(2:7), c(3:8))[, j]], c(0.025, 0.25, 0.75, 0.975), na.rm = TRUE)
      
      z_scores[which(z_scores$variable ==  c("temp", "salt")[i] & z_scores$region == region & z_scores$oe_start == j + 1), c("z_median")] <- median(z_score[, cbind(c(11:12,1:4), c(12,1:5), c(1:6), c(2:7), c(3:8))[, j]], na.rm = TRUE)
    }
    
  } # end region r
} # end variable i


write.csv(z_scores, "marine/output/NEP36-CanOE/summaryZscores_NEP36-CanOE.csv")

# Look at primary productivity