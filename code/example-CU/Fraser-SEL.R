###############################################################################
# Preliminary draft assessment of freshwater exposure factors for Fraser basin
# lake-type sockeye
# 
# Contact: Stephanie Peacock <speacock@psf.ca>
# Date: December 12, 2022
###############################################################################


library(colorRamps)
library(raster)
library(rasterVis)
library(mapproj)
library(proj4)
library(sf)
library(leaflet)
library(wesanderson)
library(ncdf4) # package for netcdf manipulation
library(zoo) # package with rollmean function

cols <- wes_palette("Darjeeling1")

###############################################################################
# Population Data
###############################################################################

#------------------------------------------------------------------------------
# List of CUs in the PSE
#------------------------------------------------------------------------------

cu_list <- read.csv("~/Documents/PSF/PopulationAnalysis/conservationunits.csv")


#------------------------------------------------------------------------------
# Timing: Spawning
#------------------------------------------------------------------------------

timing <- read.csv("data/timing/spawn-timing_CUs.csv") %>% subset(timing$region == "Fraser" & timing$SPECIES_QUALIFIED == "SEL")

timing$spawn_start <- timing$mu_spawn - 1.96*timing$sd_spawn
timing$spawn_end <- timing$mu_spawn + 1.96*timing$sd_spawn


###############################################################################
# Spatial Data Layers
###############################################################################

#------------------------------------------------------------------------------
# Conservation Unit Boundaries
#------------------------------------------------------------------------------

shp_path <- "data/spatial/CU-boundaries/Lake_Type_Sockeye_Salmon_CU_Boundary/SEL_CU_BOUNDARY_En" 

cu_boundary <- st_read(dsn = paste0(shp_path, ".shp"))
# cu_boundary <- cu_boundary[cu_boundary$CU_NAME == "Chilko-Summer Timing", ]

st_geometry_type(cu_boundary)
st_crs(cu_boundary)

# What is the spatial extent of the cu boundaries?
st_bbox(cu_boundary)

#------------------------------------------------------------------------------
# Spawning Zone of Influence
#------------------------------------------------------------------------------
zoi <- st_read(dsn = "data/spatial/ZOI/fraser-spawning-zoi/fraser_spawning_zoi_SEL.shp") %>% st_transform(crs = 4269)
#%>% subset(zoi$CUName == "Chilko_Summer") 

#------------------------------------------------------------------------------
# Lakes, rivers, and shorelines
#------------------------------------------------------------------------------

# Set bounds for clipping other files
# bounds <- c(xmin = -124.5, ymin = 50.85, xmax = -123.6, ymax = 51.7)

lakes <- readRDS("data/spatial/layers/waterbodies_250.rds")# %>% st_crop(bounds)
rivers <- readRDS("data/spatial/layers/watercourse_250.rds") #%>% st_crop(bounds)


bounds0 <- c(xmin = -127, ymin = 49, xmax = -116.5, ymax = 56)
lakes0 <- readRDS("data/spatial/layers/waterbodies_lowRes.rds")# %>% st_crop(bounds)
rivers0 <- readRDS("data/spatial/layers/watercourse_lowRes.rds") #%>% st_crop(bounds)
BC <- readRDS("data/spatial/layers/BC_lowRes.rds")
#------------------------------------------------------------------------------
# Spawning locations
#------------------------------------------------------------------------------

spawning_points <- st_read(dsn = "data/spatial/spawning_points_fraser/spawning_points_fraser.shp") %>% subset(spawning_points$Species == "lake sockeye") %>% st_transform(crs = 4269) %>% st_crop(bounds0)
# 
spawning_lines <- st_read(dsn = "data/spatial/spawning_lines_fraser/spawning_lines_fraser.shp") %>% subset(spawning_lines$Species == "lake sockeye") %>% st_transform(crs = 4269) %>% st_crop(bounds0)


###############################################################################
# PCIC *Preliminary* climate model output for the Fraser Basin
###############################################################################

#-----------------------------------------------------------------------------
# Load spatial grid: x 0.0625 degrees
#-----------------------------------------------------------------------------
link_dat <- "data/raw-data/PCIC/hydro_model_out/Fraser-prelim/"

# Subset to remove duplicate center points with V7 == 2
spat_grid <- read.csv(paste0(link_dat, "MOUTH.Spat.csv"), header = FALSE, col.names = c("V1", "V2", "V3", "V4", "lat", "lon", "V7", "V8", "V9")) %>% subset(V7 == 1)

# V2 is the unique identifier
# Note that in the prelim data there were repeat lat/lon points where two rivers meet
# Choose the one with the highest flow

# Make grid center points
z <- data.frame(lon = as.numeric(spat_grid$lon), lat = as.numeric(spat_grid$lat))
grid_points <- st_as_sf(z, coords = c("lon", "lat"), crs = 4269)

# Make grid polygon
n <- dim(spat_grid)[1]
y <- list(); length(y) <- n
d <- 0.0625
for(i in 1:n){
    y[[i]][[1]] <- rbind(
      c(spat_grid$lon[i] - d/2, spat_grid$lat[i] - d/2),
      c(spat_grid$lon[i] - d/2, spat_grid$lat[i] + d/2),
      c(spat_grid$lon[i] + d/2, spat_grid$lat[i] + d/2),
      c(spat_grid$lon[i] + d/2, spat_grid$lat[i] - d/2),
      c(spat_grid$lon[i] - d/2, spat_grid$lat[i] - d/2))
  }

grid_polys <- st_multipolygon(y) %>% st_sfc(crs = 4269)

plot(grid_polys, col = "#00000030", lwd = 0.5, border = NA)
#-----------------------------------------------------------------------------
# Open
#-----------------------------------------------------------------------------

z <- nc_open(paste0(link_dat, "MOUTH_CanESM2_r1i1p1_rcp85_sel-001.nc"))

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

# Extract temperature (Ts) and flow (Qs) variables
Ts <- ncvar_get(z, "Ts")
Qs <- ncvar_get(z, "Qs")

# Create time variable
date <- as.Date(c(0:56612), origin = "1945-01-01")

#-----------------------------------------------------------------------------
# Define decades for decade averages
#-----------------------------------------------------------------------------

decades <- data.frame(
  period = c("hist", "early", "mid", "late"),
  start = c(as.Date("1970-01-01"), as.Date("2010-01-01"), as.Date("2040-01-01"), as.Date("2070-01-01")),
  end = c(as.Date("2000-01-01"),as.Date("2040-01-01"),as.Date("2070-01-01"),as.Date("2100-01-01"))
  )

decades$start_num <- match(decades$start, date)
decades$end_num <- match(decades$end - 1, date)

#-----------------------------------------------------------------------------
# Look at optimal temperature ranges
#-----------------------------------------------------------------------------
topt <- read.csv("data/optimum-temps_BC2001.csv") %>% subset(topt$Species == "Sockeye" & topt$Life.stage == "Spawning")

###############################################################################
###############################################################################
# For each CU
###############################################################################
###############################################################################

# # Create decoder to link timing, zoi, and cu_boundaries
# write.csv(data.frame(cuid = zoi$cuid, CUName = zoi$CUName), file = "output/Fraser_SEL_zoi_decoder.csv")
# 
# write.csv(data.frame(
#   CU_NAME = cu_boundary$CU_NAME,
#   FULL_CU_IN = cu_boundary$FULL_CU_IN,
#   SP_QUAL = cu_boundary$SP_QUAL,
#   CU_TYPE = cu_boundary$CU_TYPE)
#   , file = "output/Fraser_SEL_cu_boundary_decoder.csv")
# 
#-----------------------------------------------------------------------------
# Select CU
# ** Need to streamline this with better decorder file **
#-----------------------------------------------------------------------------
# Chilko Summer CUID 721

# i <- 7
# print(timing$SQ_CU_NAME[i])
# 
# cuid <- zoi$cuid[zoi$CUName == "Chilko_Summer"] # 721
# 
# timing.i <- round(timing[which(timing$CU_NAME == "CHILKO-SUMMER TIMING"), c("spawn_start", "spawn_end")])
# zoi.i <- zoi[which(zoi$cuid == cuid), ]
# cu_boundary.i <- cu_boundary[which(cu_boundary$CU_NAME == "Chilko-Summer Timing"),]


SELcus <- read.csv("output/Fraser_SEL_decoder.csv")
n.cu <- nrow(SELcus)

# Create array to hold values
# exposure <- array(NA, dim = c(2, n.cu, 4, 2), dimnames = list(c("stream.temperature", "low.flow"), SELcus$CU_NAME_cu_boundary, decades$period, c("absolute", "percent")))

stream.temp <- list(
  NA, 
  dim = c(n.cu, 4, 2), 
  dimnames = list(
    SELcus$CU_NAME_cu_boundary, 
    decades$period, 
    c("absolute", "percent")))


for(i in 1:n.cu){
  # i = 9 = Chilko ES
  timing.i <- round(timing[which(timing$SQ_CU_NAME == paste0("SEL ", SELcus$CU_NAME[i])), c("spawn_start", "spawn_end")])
  
  zoi.i <- zoi[which(zoi$cuid == SELcus$cuid[i]), ]
  
  cu_boundary.i <- cu_boundary[which(cu_boundary$CU_NAME == SELcus$CU_NAME_cu_boundary[i]),]
  
#-----------------------------------------------------------------------------
# Identify grid cells that overlap CU
#-----------------------------------------------------------------------------

# Find those whose centers overlap CU or zoi
qqq <- st_intersects(grid_points2, zoi.i)
incl <- which(c(qqq) == 1)

# # Highlight those cells
# y_incl <- list()
# for(i in 1:length(incl)) y_incl[[i]] <- y[[incl[i]]]

# grid_polys_incl <- st_polygon(y_incl) %>% st_sfc(crs = 4269)

# # Extract relevant grid cells
# incl_ind <- cbind(match(grid_points[[1]][incl, 'lon'], lon), match(grid_points[[1]][incl, 'lat'], lat))



#------------------------------------------------------------------------------
# Plot CU
#------------------------------------------------------------------------------
if(3 == 2){
# Context map
# pdf(file = "code/example-CU/Chilko-summer/CU_context.pdf", width = 6, height = 5)
quartz(width = 6, height = 5)
plot(st_geometry(BC), axes= TRUE, , xlim = bounds0[c('xmin', 'xmax')], ylim = bounds0[c('ymin', 'ymax')])
plot(st_geometry(zoi.i), border = cols[3], col = paste0(cols[3], 50), lwd = 1.5, add = TRUE)
plot(st_geometry(cu_boundary.i), col = NA, lwd = 0.8, add = TRUE)
plot(st_geometry(lakes0), border = cols[5], col = paste0(cols[5], 50), add = TRUE)
plot(st_geometry(rivers0), col = cols[5], add = TRUE)

pdf(file = "code/example-CU/Chilko-summer/CU_map.pdf", width = 5, height = 5)
# quartz(width = 5, height = 5)

plot(st_geometry(zoi.i), border = cols[3], col = paste0(cols[3], 50), lwd = 1.5, axes= TRUE)
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
# plot(grid_polys, border = 1, lwd = 0.5, lty = 3, add = TRUE)
# plot(grid_polys_incl, col = "#00FF0030", add = TRUE)

i <- i + 1
plot(st_polygon(grid_polys_incl[[1]][i]), col = "#FF000030", add = TRUE)

plot(new_grid_polys, col = 1)

plot(st_multipoint(grid_points2[[1]][-incl,]), add = TRUE, cex = 0.6)
plot(st_geometry(st_multipoint(grid_points[[1]][incl,])), add = TRUE, cex = 0.6, pch = 19)
# dev.off()
}
#-----------------------------------------------------------------------------
# Extract time series for relevant grid cells
#-----------------------------------------------------------------------------

Ts.i <- rollmax(apply(Ts[incl, ], 2, mean, na.rm = TRUE), k = 7, align = "right", na.pad = TRUE)
# Used na.rm = TRUE because for some large ZOIs (e.g., Kamploops) there were NAs in Ts[incl,]
Qs.i <- rollmean(apply(Qs[incl, ], 2, mean, na.rm = TRUE), k = 7, align = "right", na.pad = TRUE)

# # Historical hydrograph
# p <- "late"
# h <- cbind(tapply(Qs.i[decades$start_num[decades$period == p]:decades$end_num[decades$period == p]], as.numeric(strftime(date[decades$start_num[decades$period == p]:decades$end_num[decades$period == p]], format = "%j")), median), tapply(Qs.i[decades$start_num[decades$period == p]:decades$end_num[decades$period == p]], as.numeric(strftime(date[decades$start_num[decades$period == p]:decades$end_num[decades$period == p]], format = "%j")), quantile, 0.05), tapply(Qs.i[decades$start_num[decades$period == p]:decades$end_num[decades$period == p]], as.numeric(strftime(date[decades$start_num[decades$period == p]:decades$end_num[decades$period == p]], format = "%j")), quantile, 0.95))
# plot(as.Date(1:366, origin = "2000-01-01"), h[, 1], "l", las = 1, ylab = "Flow (m/s)", xlab = "", xlim = as.Date(c(150, 330), origin = "2000-01-01"), ylim = c(0, max(h[, 3])))
# polygon(x = c(as.Date(1:366, origin = "2000-01-01"), rev(as.Date(1:366, origin = "2000-01-01"))), y = c(h[, 2], rev(h[, 3])), col = "#00000030", border = NA)
# 
# lines(as.Date(1:366, origin = "2000-01-01"), h[, 1], col = 2, lwd = 1.5)

# Take max each year during spawning period
yrs <- 1945:2099
Ts.max <- numeric(length(yrs))
Qs.min <- numeric(length(yrs))
Qs.max <- numeric(length(yrs))
for(j in 1:length(yrs)){
  Ts.max[j] <- max(Ts.i[which(date >= as.Date(paste(yrs[j], timing.i$spawn_start, sep = "-"), format = "%Y-%j") & date <= as.Date(paste(yrs[j], timing.i$spawn_end, sep = "-"), format = "%Y-%j"))])
  
  Qs.min[j] <- min(Qs.i[which(date >= as.Date(paste(yrs[j], timing.i$spawn_start, sep = "-"), format = "%Y-%j") & date <= as.Date(paste(yrs[j], timing.i$spawn_end, sep = "-"), format = "%Y-%j"))])
  
  Qs.max[j] <- max(Qs.i[which(date >= as.Date(paste(yrs[j], timing.i$spawn_start, sep = "-"), format = "%Y-%j") & date <= as.Date(paste(yrs[j], timing.i$spawn_end, sep = "-"), format = "%Y-%j"))])
}

# Calculate averages over each period
decades$Ts <- NA
for(j in 1:4){
  decades$Ts[j] <- mean(Ts.max[yrs %in% c(c(1970, 2010, 2040, 2070)[j]:c(1999, 2039, 2069, 2099)[j])])
}

exposure[1, i, , 1] <- decades$Ts
exposure[1, i, , 2] <- (decades$Ts - decades$Ts[1])/decades$Ts[1]*100

decades$Qs.low <- NA
for(j in 1:4){
  decades$Qs.low[j] <- mean(Qs.min[yrs %in% c(c(1970, 2010, 2040, 2070)[j]:c(1999, 2039, 2069, 2099)[j])])
}

# Calculate % change for each period from historical baseline
exposure[2, i, , 1] <- decades$Qs.low
exposure[2, i, , 2] <- (decades$Qs.low - decades$Qs.low[1])/decades$Qs.low[1]*100

} # end all CUs

###############################################################################
# Single CU plots
###############################################################################

#-----------------------------------------------------------------------------
# Plots: temperature
#-----------------------------------------------------------------------------
if(3 == 2){ # Nonsense if statement to block plots if sourcing script

quartz(width = 12, height = 4)

plot(date[7:length(date)], Ts.i, "l", xaxs = "i", las = 1, ylab = "Stream temperature (˚C)", xlab = "")
abline(v = decades[decades$period == "mid", c("start", 'end')], col = 2)

lines(as.Date(paste(yrs, round(mean(as.numeric(timing.i))), sep = "-"), format = "%Y-%j"), Ts.max, col = "#64837B", lwd  = 2)

segments(x0 = decades$start, x1 = decades$end, y0 = decades$Ts, y1 = decades$Ts, lwd = 4, col = "#985A3F")

quartz(width = 6, height = 4)

plot(date[decades$start_num[decades$period == "mid"]:decades$end_num[decades$period == "mid"]], Ts.i[decades$start_num[decades$period == "mid"]:decades$end_num[decades$period == "mid"]], "l", xaxs = "i", las = 1, col = grey(0.7), xlab = "", ylab = "Stream temperature (˚C)")

lines(date[c(decades$start_num[decades$period == "mid"]:decades$end_num[decades$period == "mid"])+7], Ts.i2[c(decades$start_num[decades$period == "mid"]:decades$end_num[decades$period == "mid"])+7])

quartz(width = 4, height = 4)
plot(date[7:length(date)], Ts.i, "l", xaxs = "i", las = 1, ylab = "Stream temperature (˚C)", xlab = "", xlim = date[c(decades$start_num[decades$period == "mid"]+c(0,365))])
abline(v = date[as.numeric(decades$start_num[decades$period == "mid"] + timing.i)], col = 2)

#-----------------------------------------------------------------------------
# Plots: low flow
#-----------------------------------------------------------------------------
quartz(width = 12, height = 4)

plot(date, Qs.i, "l", xaxs = "i", las = 1, ylab = "Stream flow (cm/s)", xlab = "", log = "y")
abline(v = decades[decades$period == "mid", c("start", 'end')], col = 2)

lines(as.Date(paste(yrs, round(mean(as.numeric(timing.i))), sep = "-"), format = "%Y-%j"), Qs.min, col = "#64837B", lwd  = 2)

segments(x0 = decades$start, x1 = decades$end, y0 = decades$Qs, y1 = decades$Qs, lwd = 4, col = "#985A3F")

# CLose up of mid-century
quartz(width = 6, height = 4)

plot(date[decades$start_num[decades$period == "mid"]:decades$end_num[decades$period == "mid"]], Qs.i[decades$start_num[decades$period == "mid"]:decades$end_num[decades$period == "mid"]], "l", xaxs = "i", las = 1, xlab = "", ylab = "Stream flow (cm/s)", log = "y")

#-----------------------------------------------------------------------------
# Plots: high flow
#-----------------------------------------------------------------------------
quartz(width = 12, height = 4)

plot(date, Qs.i, "l", xaxs = "i", las = 1, ylab = "Stream flow (cm/s)", xlab = "", log = "y", col = grey(0.6))
abline(v = decades[decades$period == "mid", c("start", 'end')], col = 2)

lines(as.Date(paste(yrs, round(mean(as.numeric(timing.i))), sep = "-"), format = "%Y-%j"), Qs.min, col = "#64837B50", lwd  = 2)
lines(as.Date(paste(yrs, round(mean(as.numeric(timing.i))), sep = "-"), format = "%Y-%j"), Qs.max, col = "#64837B", lwd  = 2)

segments(x0 = decades$start, x1 = decades$end, y0 = decades$Qs, y1 = decades$Qs, lwd = 4, col = "#985A3F")

# CLose up of mid-century
quartz(width = 6, height = 4)

plot(date[decades$start_num[decades$period == "mid"]:decades$end_num[decades$period == "mid"]], Qs.i[decades$start_num[decades$period == "mid"]:decades$end_num[decades$period == "mid"]], "l", xaxs = "i", las = 1, xlab = "", ylab = "Stream flow (cm/s)", log = "y")

}


###############################################################################
# Across CU plots
###############################################################################

# Temperature % change
CUcol <- PNWColors::pnw_palette("Bay", n = n.cu)

o <- order(exposure[1, , 1, 1])

quartz(width = 13, height = 8)
par(mar = c(25,4,3,2))
bp <- barplot(t(exposure[1, o, , 1]), col = rep(CUcol, each = 4), las = 2, ylim = c(0, 35), beside = TRUE)

# Relative
par(mar = c(25,4,3,2))
bp <- barplot(t(exposure[1, o, 2:4, 2]), col = rep(CUcol, each = 3), las = 2, beside = TRUE)

# Absolute change
par(mar = c(25,4,3,2))
bp <- barplot(t(exposure[1, o, 2:4, 1] - exposure[1, o, 1, 1]), col = rep(CUcol, each = 3), las = 2, beside = TRUE)


# On a map
cu_boundary2 <- cu_boundary[match(SELcus$CU_NAME_cu_boundary, cu_boundary$CU_NAME), ]

plot(st_geometry(cu_boundary2), col = CUcol[o], border = NA)
plot(st_geometry(BC), add = TRUE, col = NA)

# Low flow
quartz(width = 13, height = 8)
par(mar = c(22,4,3,2))
bp <- barplot(t(exposure[2, o, , 1]), col = rep(CUcol, each = 4), las = 2, ylim = c(0, 35), beside = TRUE, ylab = "Low flow (m/s)")

# Relative
par(mar = c(25,4,3,2))
bp <- barplot(t(exposure[2, o, 2:4, 2]), col = rep(CUcol, each = 3), las = 2, beside = TRUE)


###############################################################################
# Compare air temp to stream temp
###############################################################################

#-----------------------------------------------------------------------------
# Open max air temp
#-----------------------------------------------------------------------------

link_dat2 <- "data/raw-data/PCIC/downscaled_gcms/"
z2 <- nc_open(paste0(link_dat2, "tasmax_day_BCCAQv2+ANUSPLIN300_CanESM2_historical+rcp85_r1i1p1_19500101-21001231.nc.nc"))

# Explore dataset. What are the variables?
print(paste("File",z2$filename,"contains",z2$nvars,"variables"))
for( i in 1:z2$nvars ) {
  v <- z2$var[[i]]
  print(paste("Here is information on variable number",i))
  print(paste("   Name: ",v$name))
  print(paste("   Units:",v$units))
  print(paste("   Missing value:",v$missval))
  print(paste("   # dimensions :",v$ndims))
  print(paste("   Variable size:",v$varsize))
}

# Extract temperature (Ts) and flow (Qs) variables
lat <- ncvar_get(z2, "lat")
lon <- ncvar_get(z2, "lon")
timee <- ncvar_get(z2, "time")
length(lat)
length(lon)

# Only extract grid points that are in Fraser basin
lat.ind <- which(lat >= min(z[, 2]) & lat <= max(z[, 2]))
lon.ind <- which(lon >= min(z[, 1]) & lon <= max(z[, 1]))

tasmax <- ncvar_get(nc = z2, varid = "tasmax", start = c(min(lon.ind), min(lat.ind), 1), count = c(length(lon.ind), length(lat.ind), length(timee)))


grid_points3 <- st_as_sf(
  data.frame(
    lon = rep(lon[lon.ind], each = length(lat.ind)), 
    lat = rep(lat[lat.ind], length(lon.ind)),
    lon.ind = rep(1:length(lon.ind), each = length(lat.ind)),
    lat.ind = rep(1:length(lat.ind), length(lon.ind))
  ), 
  coords = c("lon", "lat"), crs = 4269)



# Create time variable
date2 <- as.Date(c(0:55115), origin = "1950-01-01")

#-----------------------------------------------------------------------------
# Look at relationship between change in air and stream temp for each CU
#-----------------------------------------------------------------------------

# Create array to hold values
air.exposure <- array(NA, dim = c(1, n.cu, 4, 2), dimnames = list(c("air.temperature"), SELcus$CU_NAME_cu_boundary, decades$period, c("absolute", "percent")))

for(i in 1:n.cu){
  
  timing.i <- round(timing[which(timing$CU_NAME == SELcus$CU_NAME[i]), c("spawn_start", "spawn_end")])
  
  zoi.i <- zoi[which(zoi$cuid == SELcus$cuid[i]), ]
  
  cu_boundary.i <- cu_boundary[which(cu_boundary$CU_NAME == SELcus$CU_NAME_cu_boundary[i]),]
  
  # Find those whose centers overlap CU or zoi
  qqq <- st_intersects(grid_points3, zoi.i)
  incl <- which(c(qqq) == 1)
  
  # Extract time series for relevant grid cells
  tasmax.spat <- array(NA, dim = c(length(incl), dim(tasmax)[3]))
  for(j in 1:length(incl)){
    tasmax.spat[j, ] <- tasmax[grid_points3$lon.ind[incl[j]], grid_points3$lat.ind[incl[j]], ]
  }
  
  tasmax.i <- rollmax(apply(tasmax.spat, 2, mean, na.rm = TRUE), k = 7, align = "right", na.pad = TRUE)
  # Used na.rm = TRUE because for some large ZOIs (e.g., Kamploops) there were NAs in Ts[incl,]
  
  # Take max each year during spawning period
  yrs <- 1950:2100
  tasmax.max <- numeric(length(yrs))
  for(j in 1:length(yrs)){
    tasmax.max[j] <- max(tasmax.i[which(date2 >= as.Date(paste(yrs[j], timing.i$spawn_start, sep = "-"), format = "%Y-%j") & date2 <= as.Date(paste(yrs[j], timing.i$spawn_end, sep = "-"), format = "%Y-%j"))])
  }
  
  # Calculate averages over each period
  decades$tasmax <- NA
  for(j in 1:4){
    decades$tasmax[j] <- mean(tasmax.max[yrs %in% c(c(1970, 2010, 2040, 2070)[j]:c(1999, 2039, 2069, 2099)[j])])
  }
  
  air.exposure[1, i, , 1] <- decades$tasmax
  air.exposure[1, i, , 2] <- (decades$tasmax - decades$tasmax[1])/decades$tasmax[1]*100

  }
  
#-----------------------------------------------------------------------------
# Look at relationship between air and stream temperature exposure
#-----------------------------------------------------------------------------

j <- 3 # mid

CUcol <- PNWColors::pnw_palette("Bay", n = n.cu)

quartz(width = 8, height = 3)
par(mfrow = c(1,3), mar = c(2,3,1,0), oma = c(3,3,1,1))
for(j in 2:4){
  plot(exposure[1, , j, 2], air.exposure[1, , j, 2], xlab = "", ylab = "", las = 1, pch = 21, col = CUcol, bg = paste0(CUcol, "50"))
  abline(a = 0, b = 1, lty = 3)
  abline(h = 0)
}
mtext(side = 1, "% change in stream temperature", line = 1, outer = TRUE)
mtext(side = 2, "% change in air temperature", line = 1, outer = TRUE)

plot(1, 1, "n", xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
legend("center", pch = 21, col = CUcol, pt.bg = paste0(CUcol, "50"), SELcus$CU_NAME_cu_boundary, bty = "n", xpd = NA)
legend("center", fill = CUcol, SELcus$CU_NAME_cu_boundary[o], bty = "n", xpd = NA, border = NA)

#-----------------------------------------------------------------------------
# Plot timing
#-----------------------------------------------------------------------------

x <- c(150:330)
xDate <- as.Date(x, origin = "2001-01-01")

y <- array(NA, dim = c(length(x), nrow(timing)))
for(i in 1:nrow(timing)){
  y[, i] <- dnorm(x, mean = timing$mu_spawn[i], sd = timing$sd_spawn[i])
}

plot(xDate, y[, 1], "n", ylim = range(y), yaxt = "n", ylab = "", xlab = "Spawn timing", xaxs = "i", bty = "n")
for(i in 1:n.cu){
  lines(xDate, y[, which(timing$CU_NAME == SELcus$CU_NAME[[o[i]]])], col = CUcol[i], lwd = 1.5)
  # lines(xDate, y[, which(timing$CU_NAME == SELcus$CU_NAME[[o[i]]])], col = paste0(CUcol[i], 50), lwd = 1.5)
}
for(i in c(5,7)){
  lines(xDate, y[, which(timing$CU_NAME == SELcus$CU_NAME[[o[i]]])], col = CUcol[i], lwd = 3)
}
  
# Start and end dates
plot(range(xDate), c(1,n.cu), "n", yaxt = "n", ylab = "", xlab = "Spawn timing", xaxs = "i", bty = "n")
for(i in 1:n.cu){
  timing.i <- round(timing[which(timing$CU_NAME == SELcus$CU_NAME[o][i]), c("spawn_start", "spawn_end")])
  segments(
    x0 = as.Date(as.numeric(timing.i[1]), origin = "2001-01-01"), 
    x1 = as.Date(as.numeric(timing.i[2]), origin = "2001-01-01"),
    y0 = i, y1 = i, col = CUcol[i], lwd = 7)
  # lines(xDate, y[, which(timing$CU_NAME == SELcus$CU_NAME[[o[i]]])], col = paste0(CUcol[i], 50), lwd = 1.5)
}

