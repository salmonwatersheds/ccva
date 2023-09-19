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

###############################################################################
# Load background spatial layers
###############################################################################

# Set broad bounds for clipping broad spatial layers
bounds <- c(xmin = -124.5, ymin = 50.85, xmax = -123.6, ymax = 51.7)

#------------------------------------------------------------------------------
# Lakes, rivers, and shorelines
#------------------------------------------------------------------------------

bounds0 <- c(xmin = -140, ymin = 45, xmax = -135, ymax = 63)
lakes0 <- readRDS("data/spatial/layers/waterbodies_lowRes.rds")# %>% st_crop(bounds)
rivers0 <- readRDS("data/spatial/layers/watercourse_lowRes.rds") #%>% st_crop(bounds)
BC <- readRDS("data/spatial/layers/BC_lowRes.rds")

shoreline <- st_read(dsn = "data/spatial/layers/GSHHS_i_L1.shp") %>% 
  st_transform(crs = 4269)

# buffer shoreline 400 km
shoreline400 <- st_buffer(shoreline, dist = 6.3) # Distance is in arc_degree

# Get the error: 
#   Error in wk_handle.wk_wkb(wkb, s2_geography_writer(oriented = oriented,  : 
#                                                        Loop 0 is not valid: Edge 1 crosses edge 18
# see here for solution: https://github.com/r-spatial/sf/issues/1902
sf_use_s2()
sf_use_s2(TRUE)

# Try rnaturalearth package coastline
library(rnaturalearth)
coastline <- ne_coastline(returnclass = "sf") %>%
  st_combine() %>%
  st_transform(robin)
plot(st_geometry(coastline), xlim = bounds0[c(1,3)], ylim = bounds0[c(2,4)], col = grey(0.8))

plus400 <- st_buffer(x = coastline, dist = buffer_in_m) #%>% st_crop(bounds0)

# https://gis.stackexchange.com/questions/373691/buffer-coastlines
robin <- CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

ggplot() +
  geom_sf(data = coastalWaters, fill = "lightblue", col = "transparent") +
  geom_sf(data = ROI)

#Please check the math on this part.
#KM/earth circumference * degrees in circle
buffer_in_km <- 400
buffer_in_m <- 400*1000
buffer_as_arc_degrees<- buffer_in_km/40075*360

coast <- ne_countries(returnclass = 'sf') %>%
  st_combine() %>%
  st_transform(robin)

plus400 <- coast %>%
  st_buffer(buffer_in_m)

ggplot() +
  geom_sf(data = plus400, fill = "lightblue", col = "transparent") +
  geom_sf(data = shoreline) +
  coord_sf(xlim = c(-140, -110), ylim = c(40, 60))
#------------------------------------------------------------------------------
# Adult migration lines: Fraser SEL
#------------------------------------------------------------------------------
mig_paths <- st_read(dsn = "data/spatial/fw-migration/spawn_timing_migration_paths_wCU_fraser.shp") %>% st_transform(crs = 4269)

plot(st_geometry(plus400), col = paste0(cols[5], 50), border = NA,  xlim = bounds0[c(1,3)], ylim = bounds0[c(2,4)])
plot(st_geometry(shoreline), col = grey(0.8),  xlim = bounds0[c(1,3)], ylim = bounds0[c(2,4)])
plot(st_geometry(plus400), col = paste0(cols[5], 50))

plot(st_geometry(mig_paths[which(mig_paths$cuid == 751), ]), add = TRUE, col = 2)
