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

# All NOAA grid points
NOAA_grid <- read.csv("marine/data/raw-data/NOAA/NOAA-grid-points.csv") #%>%
  # subset(hasData == 1)
# grid_points <- st_as_sf(NOAA_grid, coords = c("lon", "lat"), crs = 4269)

grid_points <- st_as_sf(NOAA_grid, coords = c("lon", "lat"), crs = 4269) %>%
  st_transform(crs = 3832)
# Use projected coordinate system to help with Pacific crossing 180 degrees

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
# Map Fraser early marine zone
###############################################################################

oe <- data.frame(lat = 49.114337, lon = -123.192161) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4269) %>%
  st_transform(crs = 3832)

fraser_1000 <- st_buffer(fraser_oe, dist = 1000*1000)
fraser_750 <- st_buffer(fraser_oe, dist = 750*1000)
fraser_300 <- st_buffer(fraser_oe, dist = 300*1000)

# Define the point where southern and northern circles should meet
breakpt <- data.frame(lat = 47.9, lon = -124.268787) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4269) %>%
  st_transform(crs = 3832) %>%
  st_bbox()

Ncrop <- st_bbox(fraser_1000)
Ncrop['ymin'] <- breakpt['ymin']
Scrop <- st_bbox(fraser_1000)
Scrop['ymax'] <- breakpt['ymax']

fraser_750N <- st_crop(fraser_750, y = Ncrop)
fraser_300S <-  st_crop(fraser_300, y = Scrop)
coastBreak <- st_union(fraser_750N, fraser_300S)

# Take difference of coastBreak and continental shelf
interior <- st_difference(x = coastBreak, y = cont)

# Take intersection of difference and pacific to cutoff land
em_zone <- st_intersection(interior, pacific)


plot(st_geometry(fraser_750), lty = 2)
plot(st_geometry(fraser_300), add = TRUE, lty = 2)
# plot(st_geometry(coastBreak), lwd = 1.5, add = TRUE)
plot(st_geometry(pacific), add = TRUE, col = paste0(cols[2], 40), border = NA)
plot(st_geometry(cont), add = TRUE, border = cols[2])
plot(st_geometry(em_zone), add = TRUE, col = paste0(cols[1], 80), border = NA)
plot(st_geometry(oe), pch = 19, add = TRUE)


# Find intersection of Pacific and fraser_300S
###############################################################################
# Spatial join of pacific shapefile to assign NOAA grid points to ocean area
###############################################################################
)

###############################################################################
# Spatial join of pacific shapefile to assign NOAA grid points to ocean area
###############################################################################

oceanRegions <- sort(unique(pacific$name))
oceanRegions_marineRearing <- c("Bering Sea", "Gulf of Alaska" , "North Pacific Ocean", "Sea of Okhotsk")

grid_joined <- st_join(grid_points, pacific)  %>%
  subset(!is.na(name)) %>% # Remove points not assigned to a region
  st_transform(crs = 4269)

# Append the ocean region to the NOAA grid data frame
NOAA_grid$region <- grid_joined$name[match(NOAA_grid$id, grid_joined$id.x)]

# Plot
plot(NOAA_grid$lon, NOAA_grid$lat, pch= 19, cex = 0.3, col = c(cols[3], grey(0.8), cols[2], rep(grey(0.8), 2), cols[5], grey(0.8), cols[4], grey(0.8), cols[1], grey(0.8))[as.numeric(factor(NOAA_grid$region, levels = names(numPoints)))])
  
write.csv(NOAA_grid, file = "marine/data/raw-data/NOAA/NOAA-grid-points_wRegion.csv")

# Plot number of grid points in each ocean area
numPoints <- tapply(grid_joined$name, grid_joined$name, length)
barplot(numPoints, col = c(cols[3], grey(0.8), cols[2], rep(grey(0.8), 2), cols[5], grey(0.8), cols[4], grey(0.8), cols[1], grey(0.8)), las = 2, ylim = c(0, 600))
plot(st_geometry(grid_joined), col = c(cols[3], grey(0.8), cols[2], rep(grey(0.8), 2), cols[5], grey(0.8), cols[4], grey(0.8), cols[1], grey(0.8))[as.numeric(factor(grid_joined$name, levels = names(numPoints)))], pch = 19, cex = 0.3)
plot(st_geometry(pacific), border = 1, col = NA, add = TRUE)

pacific_incl <- subset(pacific, pacific$name %in% c("Gulf of Alaska", "Bering Sea", "North Pacific Ocean")) %>%
  # st_transform(crs = 4269) %>%
  # st_crop(c(xmin = -180, xmax = 180, ymin = 35, ymax = 67)) %>%
  # st_transform(crs = 3832) %>%
  summarise()

# Subset grid points to North Pacific
intrscts <- st_intersects(grid_points, pacific_incl, sparse = FALSE)
incl <- which(apply(intrscts, 1, sum) == TRUE)

# Only include points north of 35 
box <- st_bbox(grid_points[incl,])
lower35 <- st_as_sf(data.frame(id= 1, lat = 35, lon = 120), coords = c("lon", "lat"), crs = 4269) %>%
  st_transform(crs = 3832)
box['ymin'] <- 4139373
# box['xmin'] <- 1.2e7
grid_points_incl <- grid_points[incl, ] %>%
  st_crop(box)

NOAA_grid_incl <- NOAA_grid[NOAA_grid$id %in% grid_points_incl$id, ]
plot(NOAA_grid_incl$lon, NOAA_grid_incl$lat, pch= 19, cex = 0.3, col = c(cols[3], grey(0.8), cols[2], rep(grey(0.8), 2), cols[5], grey(0.8), cols[4], grey(0.8), cols[1], grey(0.8))[as.numeric(factor(NOAA_grid_incl$region, levels = names(numPoints)))])

write.csv(NOAA_grid_incl, file = "marine/data/raw-data/NOAA/NOAA-grid-points_wRegion_incl.csv")


# Plotting
areaNames <- c("The Coastal Waters of Southeast Alaska and British Columbia", "Gulf of Alaska", "Bering Sea", "Sea of Okhotsk", "North Pacific Ocean")

plot(st_geometry(pacific), border = NA, col = grey(0.8), xaxs = "i", yaxs = "i", ylim = box[c('ymin', 'ymax')], xlim = c(box['xmin'], 1e7), bty = "o")
for(i in c(2,3,5)){
  # plot(st_geometry(pacific[pacific$name == areaNames[i], ]), col = cols[i], border = NA, add = TRUE)
  plot(st_geometry(pacific[pacific$name == areaNames[i], ]), col = "#FF000030", border = "#FF0000", add = TRUE)
  
  
}
legend("topleft", fill = cols[c(1:length(areaNames))], legend = c("Coastal Waters of Southeast AK + BC", "Gulf of Alaska", "Bering Sea", "Sea of Okhotsk", "North Pacific Ocean"), bty = "n", border = NA)

plot(st_geometry(grid_points), pch = 1, cex = 0.6, add = TRUE, lwd = 0.5)
plot(st_geometry(grid_points_incl), pch = 19, cex = 0.3, add = TRUE)

###############################################################################
# Ocean type chinook
###############################################################################

coast <- subset(pacific, name == "The Coastal Waters of Southeast Alaska and British Columbia")
# coast_b <- st_buffer(coast, dist = 400*1000)
coast_b <- st_buffer(coast, dist = 1000*1000)

intrscts_coast <- st_intersects(grid_points, coast_b, sparse = FALSE)
incl_coast <- which(intrscts_coast == TRUE)

grid_points_incl_coast <- grid_points[incl_coast, ]

plot(st_geometry(pacific), border = NA, col = grey(0.8), xaxs = "i", yaxs = "i", ylim = box[c('ymin', 'ymax')], xlim = c(box['xmin'], 1e7))

plot(st_geometry(coast), col = cols[1], border = NA, add = TRUE)
plot(st_geometry(coast_b), col = NA, border = cols[1], add = TRUE)
plot(st_geometry(grid_points_incl_coast[which(grid_points_incl_coast$hasData == 1), ]), pch = 19, cex = 0.3, add = TRUE)

plot(st_geometry(pacific_incl), col = cols[1], border = NA, add = TRUE)
plot(st_geometry(subset(pacific, name %in% c("North Pacific Ocean", "Gulf of Alaska"))), col = cols[1], border = NA, add = TRUE)
plot(st_geometry(subset(pacific, name %in% c("Bering Sea"))), col = cols[1], border = NA, add = TRUE)

NOAA_grid_incl_coast <- NOAA_grid[NOAA_grid$id %in% grid_points_incl_coast$id, ]
write.csv(NOAA_grid_incl_coast, file = "marine/data/raw-data/NOAA/NOAA-grid-points_wRegion_incl_coast1000.csv")
