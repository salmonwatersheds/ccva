###############################################################################
#
# 1a_subset-NOAA-grid.R
#
# This code reads in a number of spatial datasets related to determine which
# grid points for the coarse GCM model output provided by NOAA should be included
# for the Fraser early marine and marine rearing stages.
#
# Contact: Steph Peacock (speacock@psf.ca)
# Date: Jan 31, 2024
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

###############################################################################
# Load background spatial layers
###############################################################################

# All NOAA grid points
grid_points0 <- read.csv("marine/data/raw-data/NOAA/NOAA-grid-points.csv") #%>%
  # subset(hasData == 1)
# grid_points <- st_as_sf(NOAA_grid, coords = c("lon", "lat"), crs = 4269)

grid_points <- st_as_sf(grid_points0, coords = c("lon", "lat"), crs = 4269) %>%
  st_transform(crs = 3832)
# Use projected coordinate system to help with Pacific crossing 180 degrees

n <- length(grid_points$id)
d <- 1
grid_polys <- st_as_sf(data.frame(
  id = rep(grid_points0$id, each = 5),
  rep(c("SW0", "NW", "NE", "SE", "SW1"), n),
  lon = c(rep(grid_points0$lon, each = 5) + rep(c(- d/2, -d/2, d/2,  d/2, -d/2), n)),
  lat = c(rep(grid_points0$lat, each = 5) + rep(c(- d/2,  d/2, d/2, -d/2, -d/2), n))
), coords = c("lon", "lat"), crs = 4269) %>% 
  group_by(id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")  %>%
  st_transform(crs = 3832)
#------------------------------------------------------------------------------
# Read in shapefile of North Pacific regions
#------------------------------------------------------------------------------
pacific <- st_read(dsn = "marine/data/spatial/iho/iho.shp") %>% 
  # st_transform(crs = 4269)
  st_transform(crs = 3832)
sf_use_s2(FALSE)

#------------------------------------------------------------------------------
# Read in shapefile of Continental margin
#------------------------------------------------------------------------------

cont <- st_read(dsn = "marine/data/spatial/Global Margin/ContinentalMargins.shp") %>% 
  # st_transform(crs = 4269)
  st_transform(crs = 3832)

###############################################################################
# Early marine stage: radius around ocean entry and 100 km off ocean shelf
###############################################################################

# Define ocean entry (oe) point
fraser_oe <- data.frame(lat = 49.114337, lon = -123.192161) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4269) %>%
  st_transform(crs = 3832)

fraser_2000 <- st_buffer(fraser_oe, dist = 2000*1000)
fraser_1000 <- st_buffer(fraser_oe, dist = 1000*1000)
fraser_750 <- st_buffer(fraser_oe, dist = 750*1000)
fraser_300 <- st_buffer(fraser_oe, dist = 300*1000)

# crop cont to fraser_2000
cont_2000 <- st_crop(cont, fraser_2000)
cont_2000_buff <- st_buffer(cont_2000, dist = -100*1000)

# Define the point where southern and northern circles should meet
breakpt <- data.frame(lat = 5990419, lon = 9536413) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 3832) %>%
  st_bbox()


Ncrop <- st_bbox(fraser_1000)
Ncrop['ymin'] <- breakpt['ymin']
Scrop <- st_bbox(fraser_1000)
Scrop['ymax'] <- breakpt['ymax']

fraser_750N <- st_crop(fraser_750, y = Ncrop)
fraser_300S <-  st_crop(fraser_300, y = Scrop)
coastBreak <- st_union(fraser_750N, fraser_300S)
# # addCorner <- data.frame(ID = c(1:4), lat = c(6755951, 6755951,6285305,6285305),  lon = c(9527219, 8737202, 8737202, 9527219)) %>%
#   st_as_sf(coords = c("lon", "lat"), crs = 3832) %>%
#   group_by(ID) %>% 
#   aggregate(by = "ID") %>%
#   st_cast("POLYGON")
# 
#   plot(st_geometry(coastBreak), add = TRUE)

# Take difference of coastBreak and continental shelf
interior <- st_difference(x = coastBreak, y = cont_2000_buff)

# Take intersection of difference and pacific to cutoff land
em_zone <- st_intersection(interior, pacific)


plot(st_geometry(fraser_750), lty = 2)
plot(st_geometry(fraser_300), add = TRUE, lty = 2)
# plot(st_geometry(coastBreak), lwd = 1.5, add = TRUE)
plot(st_geometry(pacific), add = TRUE, col = paste0(cols[2], 40), border = NA)
# plot(st_geometry(cont), add = TRUE, border = cols[2], lty = 3)
plot(st_geometry(cont_2000_buff), add = TRUE, border = cols[2])
plot(st_geometry(em_zone), add = TRUE, col = paste0(cols[1], 80), border = NA)
plot(st_geometry(fraser_oe), pch = 19, add = TRUE)
plot(st_geometry(grid_points), cex = 0.3, add = TRUE, col = grey(0.7))
plot(st_geometry(grid_polys), col = NA, lwd = 0.5, add = TRUE, border = grey(0.7))
plot(st_geometry(grid_points[grid_points$hasData == 1, ]), cex = 0.3, add = TRUE)

# saveRDS(em_zone, file = "marine/output/fraser_EM_zone.rds")

#------------------------------------------------------------------------------
# determine NOAA grid cells that intersect in em_zone and have data
#------------------------------------------------------------------------------

grid_polys_in <- st_intersection(grid_polys[grid_polys$id %in% grid_points$id[grid_points$hasData == 1], ], em_zone)
id_in <- unique(grid_polys_in$id)

plot(st_geometry(grid_points[grid_points$id %in% id_in, ]), col = cols[3], cex = 1.5, pch = 19, add = TRUE)

write.csv(id_in, "marine/output/incl_NOAA_EM.csv")

###############################################################################
# Marine rearing
###############################################################################

oceanRegions <- sort(unique(pacific$name))
oceanRegions_marineRearing <- c("Bering Sea", "Gulf of Alaska" , "North Pacific Ocean")

southBound <- data.frame(lat = 35, lon = -135) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4269) %>%
  st_transform(crs = 3832)

pacific_incl <- subset(pacific, pacific$name %in% oceanRegions_marineRearing) %>%
  st_union()

# Subset grid points to North Pacific
intrscts <- st_intersects(grid_points, pacific_incl, sparse = FALSE)
incl <- which(apply(intrscts, 1, sum) == TRUE)

# Only include points north of 35 
incl_N35 <- grid_points0$id[which(grid_points0$id %in% incl & grid_points$hasData == 1 & grid_points0$lat >= 35)]

write.csv(incl_N35, file = "marine/output/incl_NOAA_marRear.csv")

