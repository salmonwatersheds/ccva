#------------------------------------------------------------------------------
# Lakes, rivers, and shorelines
#------------------------------------------------------------------------------

# lakes <- readRDS("data/spatial/layers/waterbodies_250.rds")# %>% st_crop(bounds)
# rivers <- readRDS("data/spatial/layers/watercourse_250.rds") #%>% st_crop(bounds)
# 
# bounds0 <- c(xmin = -127, ymin = 49, xmax = -116.5, ymax = 56)
lakes0 <- readRDS("data/spatial/layers/waterbodies_lowRes.rds")# %>% st_crop(bounds)
rivers0 <- readRDS("data/spatial/layers/watercourse_lowRes.rds") #%>% st_crop(bounds)
BC <- readRDS("data/spatial/layers/BC_lowRes.rds")
shoreline <- st_read(dsn = "data/spatial/layers/GSHHS_i_L1.shp")

#-----------------------------------------------------------------------------
# map life stages
#-----------------------------------------------------------------------------

# Use Sam's palette
colStages <- c(adult_migration = "#001289",
               spawning = "#0F78B6",
               eggs_alevin = "#53B3E9",
               fw_rearing = "#009C70")
z <- nc_open(paste0("data/raw-data/PCIC/hydro_model_out/fraser/", "waterTemp", "_day_dynWat-VICGL_", "CanESM2", "_rcp45_", "r1i1p1", "_19450101-20991231_fraser.nc"))

# Spatial variables
lon <- ncvar_get(z, "lon")
lat <- ncvar_get(z, "lat")

grid_points <- data.frame(
  id = c(1:(length(lon)*length(lat))),
  lon = rep(lon, length(lat)),
  lat = rep(lat, each = length(lon)),
  lon_id = rep(c(1:length(lon)), length(lat)),
  lat_id = rep(c(1:length(lat)), each = length(lon)))

n <- length(grid_points$id)
d <- 1/16

grid_polys_all <- st_as_sf(data.frame(
  id = rep(grid_points$id, each = 5),
  rep(c("SW0", "NW", "NE", "SE", "SW1"), n),
  lon = c(rep(rep(lon, length(lat)), each = 5) + rep(c(- d/2, -d/2, d/2,  d/2, -d/2), n)),
  lat = c(rep(rep(lat, each = length(lon)), each = 5) + rep(c(- d/2,  d/2, d/2, -d/2, -d/2), n))
), coords = c("lon", "lat"), crs = 4269) %>% 
  group_by(id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") 


#-----------------------------------------------------------------------------
# Select CU
#-----------------------------------------------------------------------------
i <- 19

cuid <- cu_list$CUID[i]
cu_boundary.i <- cu_boundary[which(cu_boundary$FULL_CU_IN == cu_list$Full.CU.Index[i]),]
zoi.i <- spawn_zoi[which(spawn_zoi$cuid == cuid), ]
mig_paths.i <- mig_paths[which(mig_paths$cuid == cuid), ]

#-----------------------------------------------------------------------------
# Base map of different CU-specific spatial data
#-----------------------------------------------------------------------------
# Set bounds for map as extent of migration and spawning ZOI
bounds0 <- st_bbox(mig_paths.i)[c(1,3,2,4)]
bounds1 <- st_bbox(zoi.i)[c(1,3,2,4)]
bounds2 <- apply(cbind(bounds0, bounds1), 1, range)
bounds <- c(min(bounds2[,1]), max(bounds2[,2]), min(bounds2[,3]), max(bounds2[,4]))

# Plot map
par(mar = c(2,2,2,1), bg = "white")

# Set up blank plot
plot(st_geometry(BC), border = NA, col = "white", axes = FALSE, las = 1, ylim = bounds[3:4], xlim = bounds[1:2], bty = "o")

# Add shoreline
plot(shoreline, add = TRUE, col = NA, border = 1)

# coarser lakes and rivers
plot(st_geometry(lakes0), border = grey(0.6), col = NA, add = TRUE)
plot(st_geometry(rivers0), col = grey(0.6), add = TRUE)

plot(st_geometry(grid_polys_all), border = "#00000020", add = TRUE)

plot(st_geometry(mig_paths.i), col = colStages['adult_migration'], lwd = 2, add = TRUE)
plot(st_geometry(zoi.i), border = colStages['spawning'], col = NA, lwd = 2, add = TRUE)
plot(st_geometry(cu_boundary.i), border = colStages['fw_rearing'], col = NA, lwd = 2, add = TRUE)

#-----------------------------------------------------------------------------
# Map selected grid cells
#-----------------------------------------------------------------------------
# Plot map
par(mar = c(2,2,2,1), bg = "white")

# Set up blank plot
plot(st_geometry(BC), border = NA, col = "white", axes = FALSE, las = 1, ylim = bounds[3:4], xlim = bounds[1:2], bty = "o")

# Add shoreline
plot(st_geometry(shoreline), add = TRUE, col = NA, border = 1)

# coarser lakes and rivers
plot(st_geometry(lakes0), border = grey(0.6), col = NA, add = TRUE)
plot(st_geometry(rivers0), col = grey(0.6), add = TRUE)

plot(st_geometry(grid_polys_all), border = "#00000020", add = TRUE)

#-----
# Adult migration
#-----
plot(st_geometry(mig_paths.i), col = colStages['adult_migration'], add = TRUE)
# plot(st_geometry(cu_boundary.i), border = colStages['adult_migration'], col = NA, add = TRUE)
plot(st_geometry(grid_polys_all[incl.stages[[i,1]], ]), border = NA, col = paste0(colStages['adult_migration'], 60), add = TRUE)
mtext(side = 3, adj = 0, line = -2, "   Adult freshwater migration", col = colStages['adult_migration'], font = 2, cex = 1.2)

#-----
# Spawning
#-----
plot(st_geometry(zoi.i), border = colStages['spawning'], col = NA, add = TRUE)
plot(st_geometry(grid_polys_all[incl.stages[[i,2]], ]), col = paste0(colStages['spawning'], 60), border = NA, add = TRUE)
mtext(side = 3, adj = 0, line = -2, "   Spawning", col = colStages['spawning'], font = 2, cex = 1.2)

#-----
# Incubation
#-----
plot(st_geometry(zoi.i), border = colStages['spawning'], col = NA, add = TRUE)
plot(st_geometry(grid_polys_all[incl.stages[[i,3]], ]), col = paste0(colStages['eggs_alevin'], 60), border = NA, add = TRUE)
mtext(side = 3, adj = 0, line = -2, "   Incubation", col = colStages['eggs_alevin'], font = 2, cex = 1.2)

#-----
# Freshwater rearing
#-----
plot(st_geometry(mig_paths.i), col = colStages['adult_migration'], add = TRUE)
plot(st_geometry(cu_boundary.i), border = colStages['fw_rearing'], col = NA, add = TRUE)
plot(st_geometry(grid_polys_all[incl.stages[[i,4]], ]), col = paste0(colStages['fw_rearing'], 60), border = NA, add = TRUE)
mtext(side = 3, adj = 0, line = -2, "   Freshwater rearing", col = colStages['fw_rearing'], font = 2, cex = 1.2)
