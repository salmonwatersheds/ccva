#------------------------------------------------------------------------------
# Create common list of spatial distribution: which grid cells for which stage
#------------------------------------------------------------------------------
stages <- c("eggs_alevin", "fw_rearing", "early_marine", "marine_rearing", "adult_migration", "spawning")

# Load grid cell ids for each CU and life stage
incl.stages_fw <- readRDS("data/PCIC_incl.rds")
incl.stages_fw[[1,1]]

# Load marine grid points
incl.stages_em <- read.csv("marine/output/incl_NOAA_EM.csv")
incl.stages_mr <- read.csv("marine/output/incl_NOAA_marRear.csv")

# Combine incl.stages into one

incl.stages <- list(); length(incl.stages) <- 60*6; dim(incl.stages) <- c(60, 6)
dimnames(incl.stages) <- list(dimnames(incl.stages_fw)[[1]], stages)
for(i in 1:60){
  incl.stages[[i, 1]] <- incl.stages_fw[[i, match(stages[1], dimnames(incl.stages_fw)[[2]])]]
  incl.stages[[i, 2]] <- incl.stages_fw[[i, match(stages[2], dimnames(incl.stages_fw)[[2]])]]
  incl.stages[[i, 3]] <- incl.stages_em$x
  incl.stages[[i, 4]] <- incl.stages_mr$x
  incl.stages[[i, 5]] <- incl.stages_fw[[i, match(stages[5], dimnames(incl.stages_fw)[[2]])]]
  incl.stages[[i, 6]] <- incl.stages_fw[[i, match(stages[6], dimnames(incl.stages_fw)[[2]])]]
}

saveRDS(incl.stages, file = "output/incl.stages_all.rds")

#------------------------------------------------------------------------------
# Create spatial grid cells for marine
#------------------------------------------------------------------------------

incl.stage_marAll <- sort(unique(c(incl.stages_em$x, incl.stages_mr$x)))

grid_points.mar <- read.csv("marine/data/raw-data/NOAA/NOAA-grid-points.csv") %>%
  filter(id %in% incl.stage_marAll) # 2019 x 6

# Visual check
# plot(grid_points.mar$lon[which(grid_points.mar$id %in% incl.stages_mr$x)], grid_points.mar$lat[which(grid_points.mar$id %in% incl.stages_mr$x)])

grid_points.mar1 <- st_as_sf(grid_points.mar, coords = c("lon", "lat"), crs = 4269)

plot(st_geometry(grid_points.mar1), axes = TRUE)

n <- length(grid_points.mar$id)
d <- 1
grid_polys.mar <- st_as_sf(data.frame(
  id = rep(grid_points.mar$id, each = 5),
  rep(c("SW0", "NW", "NE", "SE", "SW1"), n),
  lon = c(rep(grid_points.mar$lon, each = 5) + rep(c(- d/2, -d/2, d/2,  d/2, -d/2), n)),
  lat = c(rep(grid_points.mar$lat, each = 5) + rep(c(- d/2,  d/2, d/2, -d/2, -d/2), n))
), coords = c("lon", "lat"), crs = 4269) %>% 
  group_by(id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") # %>%
  # st_transform(crs = 3832)

plot(st_geometry(grid_polys.mar[2000:2019,]), col = NA, border = 2, add = TRUE)

saveRDS(grid_polys.mar, file = "docs/data/grid_polys_mar.rds")
