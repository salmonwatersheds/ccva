# Get province borders and project it to same CRS than raster
# From https://stackoverflow.com/questions/49521933/add-graticule-to-map-using-levelplot-r

require(colorRamps)
require(raster)
library(rasterVis)
require(mapproj)
library(raster)    
library(proj4)
library(sf)


#------------------------------------------------------------------------------
# Baselayers
#------------------------------------------------------------------------------

base_path <- "~/Documents/Mapping/lpr_000b21a_e/lpr_000b21a_e"
coastline0 <- st_read(dsn = paste0(base_path, ".shp"))

# Change coordinate system to match CU boundaries: NAD83
coastline <- st_transform(coastline0, crs = 4269)
st_crs(coastline)

names(coastline)
unique(coastline$PRENAME)
coastline <- coastline[coastline$PRENAME %in% c("British Columbia"), ]

st_geometry_type(coastline)

coastline_simplest <- st_simplify(coastline, dTolerance = 1000)
coastline_simpler <- st_simplify(coastline, dTolerance = 200)

plot(st_geometry(coastlineS), add = TRUE, col = 2)

# Rivers
rivers0 <- st_read(dsn = "~/Documents/Mapping/lhy_000d16a_e/lhy_000d16a_e.shp")
st_crs(rivers0)
rivers <- st_transform(rivers0, crs = 4269)

rivers_simplest <- st_simplify(rivers, dTolerance = 1000)
rivers_simpler <- st_simplify(rivers, dTolerance = 10)
rivers_cu <- st_intersection(rivers, cu_boundary)
rivers_simplest_cu <- st_intersection(rivers_simplest, cu_boundary)

# Lakes
lakes0 <- st_read

#------------------------------------------------------------------------------
# Conservation Unit Boundary
#------------------------------------------------------------------------------

shp_path <- "raw-data/conservation-units/Chinook/Chinook_Salmon_CU_Shape/Chinook_Salmon_CU_Boundary/CK_CU_BOUNDARY_En"

cu_boundary <- st_read(dsn = paste0(shp_path, ".shp"))

st_geometry_type(cu_boundary)
st_crs(cu_boundary)

# What is the spatial extent of the cu boundary?library
st_bbox(cu_boundary)

# Interested in "North Thompson_SU_1.3"
cu_boundary <- cu_boundary[cu_boundary$CU_NAME == "North Thompson_SU_1.3", ]

# Context plot
plot(st_geometry(coastline_simplest), axes = TRUE)
plot(st_geometry(cu_boundary), add = TRUE, col = 2)

plot(st_geometry(cu_boundary))
plot(st_geometry(rivers_simplest_cu), add = TRUE)
plot(st_geometry(rivers_cu), add = TRUE, col = 4)


# CU plot