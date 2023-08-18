# Example CU: Chilko Summer-timing lake-type sockeye

library(colorRamps)
library(raster)
library(rasterVis)
library(mapproj)
library(proj4)
library(sf)
library(leaflet)
library(wesanderson)

cols <- wes_palette("Darjeeling1")
###############################################################################
# Spatial Data Layers
###############################################################################

#------------------------------------------------------------------------------
# Conservation Unit Boundary
#------------------------------------------------------------------------------

shp_path <- "data/raw-data/conservation-units/Lake_Type_Sockeye/Lake_Type_Sockeye_Salmon_CU_Shape/Lake_Type_Sockeye_Salmon_CU_Boundary/SEL_CU_BOUNDARY_En" 

cu_boundary <- st_read(dsn = paste0(shp_path, ".shp"))
cu_boundary <- cu_boundary[cu_boundary$CU_NAME == "Chilko-Summer Timing", ]

st_geometry_type(cu_boundary)
st_crs(cu_boundary)

# What is the spatial extent of the cu boundary?library
st_bbox(cu_boundary)

#------------------------------------------------------------------------------
# Lakes, rivers, and shorelines
#------------------------------------------------------------------------------

# Set bounds for clipping other files
bounds <- c(xmin = -124.5, ymin = 50.85, xmax = -123.6, ymax = 51.7)

lakes <- readRDS("data/mapping/layers/waterbodies_250.rds") %>% st_crop(bounds)

rivers <- readRDS("data/mapping/layers/watercourse_250.rds") %>% st_crop(bounds)

#------------------------------------------------------------------------------
# Spawning Zone of Influence
#------------------------------------------------------------------------------
zoi <- st_read(dsn = "code/example-CU/Chilko-summer/fraser-spawning-zoi/fraser_spawning_zoi_SEL.shp") %>% subset(zoi$CUName == "Chilko_Summer") %>% st_transform(crs = 4269)

#------------------------------------------------------------------------------
# Spawning locations
#------------------------------------------------------------------------------

spawning_points <- st_read(dsn = "data/raw-data/conservation-units/spawning_points_fraser/spawning_points_fraser.shp") %>% subset(spawning_points$Species == "lake sockeye") %>% st_transform(crs = 4269) %>% st_crop(bounds)

spawning_lines <- st_read(dsn = "data/raw-data/conservation-units/spawning_lines_fraser/spawning_lines_fraser.shp") %>% subset(spawning_lines$Species == "lake sockeye") %>% st_transform(crs = 4269) %>% st_crop(bounds)

###############################################################################
# PCIC climate model output
###############################################################################

#-----------------------------------------------------------------------------
# Open
#-----------------------------------------------------------------------------

z <- nc_open("code/example-CU/Chilko-summer/allwsbc.ACCESS1-0_rcp85_r1i1p1.1945to2099.BASEFLOW.nc.nc")

lat <- ncvar_get(z, "lat")
lon <- ncvar_get(z, "lon")
time <- ncvar_get(z, "time")

tunits <- ncatt_get(z,"time","units")

date <- as.Date(time, origin = "1945-01-01")

bf <- ncvar_get(z, "BASEFLOW")
bf[bf == -32767] <- NA # Replace FillValues with NA

#-----------------------------------------------------------------------------
# Define decades for decade averages
#-----------------------------------------------------------------------------

decades <- data.frame(period = c("hist", "early", "mid", "late"),
                      start = c(as.Date("1970-01-01"), as.Date("2010-01-01"), as.Date("2040-01-01"), as.Date("2070-01-01")),
                      end = c(as.Date("2000-01-01"),as.Date("2040-01-01"),as.Date("2070-01-01"),as.Date("2100-01-01")))

decades$start_num <- match(decades$start, date)
decades$end_num <- match(decades$end - 1, date)

#-----------------------------------------------------------------------------
# Identify grid cells that overlap CU
#-----------------------------------------------------------------------------

# Make grid polygon
n <- length(lon)*length(lat)
d <- diff(lon)[1]/2 # lon and lat diff are the same = 0.0625 degrees
y <- list(); length(y) <- n
for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    y[[(i - 1)*(length(lat)) + j]] <- rbind(
      c(lon[i] - d, lat[j] - d),
      c(lon[i] - d, lat[j] + d),
      c(lon[i] + d, lat[j] + d),
      c(lon[i] + d, lat[j] - d),
      c(lon[i] - d, lat[j] - d))
  }
}
grid_polys <- st_polygon(y) %>% st_sfc(crs = 4269)

# Make grid center points
z <- cbind(lon = rep(lon, each = length(lat)), lat = rep(lat, length(lon)))
grid_points <- st_multipoint(z) %>% st_sfc(crs = 4269)


# Find those whose centers overlap CU or zoi
incl <- c()
for(i in 1:n){
  q <- st_intersects(st_sfc(st_point(z[i,]), crs = 4269), zoi)#cu_boundary)
  if(length(q[[1]]) == 1) incl <- c(incl, i)
}

# Highlight those cells
y_incl <- list()
for(i in 1:length(incl)) y_incl[[i]] <- y[[incl[i]]]

grid_polys_incl <- st_polygon(y_incl) %>% st_sfc(crs = 4269)

# Extract relevant grid cells
incl_ind <- cbind(match(z[incl, 'lon'], lon), match(z[incl, 'lat'], lat))

#------------------------------------------------------------------------------
# Plot CU
#------------------------------------------------------------------------------


pdf(file = "code/example-CU/Chilko-summer/CU_map4.pdf", width = 5, height = 8)
plot(st_geometry(zoi), border = cols[3], col = paste0(cols[3], 50), lwd = 1.5, axes= TRUE)
plot(st_geometry(cu_boundary), border = cols[4], col = paste0(cols[4], 50), lwd = 1.5, add = TRUE)
plot(st_geometry(lakes), border = cols[5], col = paste0(cols[5], 50), add = TRUE)
plot(st_geometry(rivers), col = cols[5], add = TRUE)
plot(st_geometry(spawning_points), col = cols[1], pch = 19, add = TRUE, cex = 0.7)
plot(st_geometry(spawning_lines), col = cols[1], add = TRUE)

# Add grid to map
plot(grid_polys, border = 1, lwd = 0.5, lty = 3, add = TRUE)
plot(grid_polys_incl, col = "#00000030", add = TRUE)
text(z[incl, 1], z[incl, 2], incl, cex = 0.6)
dev.off()

#-----------------------------------------------------------------------------
# Extract time series from one grid cell
#-----------------------------------------------------------------------------

bf2 <- bf[match(z[150, 'lon'], lon), match(z[150, 'lat'], lat), ]

pdf(file = "code/example-CU/Chilko-summer/timeseries2.pdf", width = 10, height = 4)
plot(date, bf2, "n", ylab = "Baseflow (mm)", xlab = "", las = 1, ylim = c(0, 1))
for(i in 1:length(incl)) lines(date, bf[incl_ind[i,1], incl_ind[i,2], ], col = grey(0.8), lwd = 0.8)
lines(date, bf2, col = 1, lwd = 1.2)
dev.off()

#-----------------------------------------------------------------------------
# Mean across all selected grid cells
#-----------------------------------------------------------------------------

bf_incl <- array(NA, dim = c(length(incl), length(time)))
for(i in 1:length(incl)){
  bf_incl[i, ] <- bf[incl_ind[i,1], incl_ind[i,2], ]
}

bf_avg <- cbind(mean = apply(bf_incl, 2, mean),
                lwr = apply(bf_incl, 2, quantile, 0.025),
                upr = apply(bf_incl, 2, quantile, 0.975))
                

plot(date, bf_avg[, 'mean'], "n", las = 1, ylab = "Baseflow (mm)", xlab = "", ylim = c(0, 5))
polygon(x = c(date, rev(date)), y= c(bf_avg[, 'lwr'], rev(bf_avg[, "upr"])), border = NA, col = grey(0.8))
lines(date, bf_avg[, 'mean'])


bf_decade <- c(
  hist = mean(bf_avg[decades$start_num[1]:decades$end_num[1], 'mean']),
  early = mean(bf_avg[decades$start_num[2]:decades$end_num[2], 'mean']),
  mid = mean(bf_avg[decades$start_num[3]:decades$end_num[3], 'mean']),
  late = mean(bf_avg[decades$start_num[4]:decades$end_num[4], 'mean'])
)

segments(x0 = as.Date("1970-01-01"), x1 = as.Date("2000-01-01"), y0 = bf_decade[1], y1 = bf_decade[1], col= cols[1], lwd = 2)

segments( x0 = as.Date("2010-01-01"), x1 = as.Date("2040-01-01"), y0 = bf_decade[2], y1 = bf_decade[2], col= cols[1], lwd = 2)

segments(x0 = as.Date("2040-01-01"), x1 = as.Date("2070-01-01"), y0 = bf_decade[3], y1 = bf_decade[3], col= cols[1], lwd = 2)

segments(x0 = as.Date("2070-01-01"), x1 = as.Date("2100-01-01"), y0 = bf_decade[4], y1 = bf_decade[4], col= cols[1], lwd = 2)

#-----------------------------------------------------------------------------
# Average for min and max flow
#-----------------------------------------------------------------------------

# Low flows: number of days below historical 20% MAD
mad <- apply(bf_incl[, decades$start_num[1]:decades$end_num[1]], 1, mean)

hist(0.1*mad, xlab = "10% MAD for 1970-2000", main = "Distribution among grid cells")

# Number of years in output
ny <- as.numeric(round(diff(range(date))/365))

# Define function to determine number of days x that are less than threshold y
lw <- function(x, y){
  length(which(x < y))
}

nd_below_10pmad <- array(NA, dim = c(length(incl), ny), dimnames = list(NULL, c(1945:2099)))
for(i in 1:length(incl)){ # for each grid cell
  nd_below_10pmad[i, ] <- tapply(bf_incl[i, ], as.numeric(strftime(date, format = "%Y")), lw, 0.1*mad[i])
}

ndJJA_below_10pmad <- array(NA, dim = c(length(incl), ny), dimnames = list(NULL, c(1945:2099)))
for(i in 1:length(incl)){ # for each grid cell
  ndJJA_below_10pmad[i, ] <- tapply(
    bf_incl[i, as.numeric(strftime(date, format = "%m")) %in% c(6:8)], 
    as.numeric(strftime(date[as.numeric(strftime(date, format = "%m")) %in% c(6:8)], format = "%Y")), lw, 0.1*mad[i])
}

y <- nd_below_10pmad
plot(c(1945:2099), y[1, ], "n", ylim = range(y), ylab = "Number of days < 10%MAD", xlab = "")
for(i in 1:length(incl)) lines(c(1945:2099), y[i, ], col = grey(0.8), lwd = 0.8)
# lines(c(1945:2099), y[which(incl == 150), ])
# lines(c(1945:2099), y[which(incl == 181), ])

# Summarize by decade
avg_nd_below_10pmad <- c(
  hist = mean(nd_below_10pmad[, which(c(1945:2099) %in% c(1970:1999))]),
  early = mean(nd_below_10pmad[, which(c(1945:2099) %in% c(2010:2039))]),
  mid = mean(nd_below_10pmad[, which(c(1945:2099) %in% c(2040:2069))]),
  late = mean(nd_below_10pmad[, which(c(1945:2099) %in% c(2070:2099))]))
                             
for(j in 1:4) segments(x0 = c(1970, 2010, 2040, 2070)[j], x1 = c(1999, 2039, 2069, 2099)[j], y0 = avg_nd_below_10pmad[j], y1 = avg_nd_below_10pmad[j], lwd = 2, col = cols[1])


hist(nd_below_10pmad[, which(c(1945:2099) %in% c(2040:2069))], main = "Distribution of #days \nacross grid cells and years 2040-2069", breaks = seq(0, max(nd_below_10pmad), 7), xlab = "# days", xlim = c(0, 150))
abline(v = avg_nd_below_10pmad['mid'], lwd = 2, col = cols[1])
hist(nd_below_10pmad[which(incl == 181), which(c(1945:2099) %in% c(2040:2069))], col = cols[1], border = NA, add = TRUE, breaks = seq(0, max(nd_below_10pmad), 7), freq = FALSE)

yy <- apply(nd_below_10pmad[, which(c(1945:2099) %in% c(2040:2069))], 1, mean)
xx <- apply(nd_below_10pmad[, which(c(1945:2099) %in% c(2040:2069))], 2, mean)

hist(yy, main = "Distribution of mean #days in mid-century \nacross grid cells", breaks = seq(0, max(nd_below_10pmad), 7), xlab = "# days", xlim = c(0, 150))

hist(xx, main = "Distribution of mean #days across grid cells \n among years", breaks = seq(0, max(nd_below_10pmad), 7), xlab = "# days", xlim = c(0, 150))
