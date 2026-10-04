###############################################################################
# Code to produce Figure 1 for paper: map of Fraser River basin with points
# at centroid of each CU, shaded according to overall exposure.
#
# December 23, 2025
###############################################################################
library(sf)

load("output/mapping_2025-12-03.RData")

# Filter CU boundaries to 
cu_boundary <- cu_boundary %>%
  filter(cu_boundary$CUID %in% unique(exposure$cuid))

# Order by cu_order
cu_boundary <- cu_boundary[match(cu_order, cu_boundary$CUID), ]

cu_union <- st_union(cu_boundary)
# sf_use_s2(FALSE)

# Clip the lines to the polygon
clipped_rivers <- st_intersection(rivers0, cu_union)
clipped_lakes <- st_intersection(lakes0, cu_union)

# Move centroid points to avoid overlap
cu_centroid <- st_centroid(cu_boundary)

# Function 
# Enforce minimum separation between point features
enforce_min_distance <- function(points_sf, min_dist, max_iter = 100, step = 0.1) {
  stopifnot(inherits(points_sf, "sf"))
  geom_type <- unique(sf::st_geometry_type(points_sf))
  if (!all(geom_type == "POINT")) stop("Input must be POINT geometries")
  
  coords <- sf::st_coordinates(cu_centroid)
  coords
  n <- nrow(coords)
  
  for (iter in seq_len(max_iter)) {
    moved <- FALSE
    
    # Pairwise distance matrix
    dmat <- as.matrix(dist(coords))
    
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        d <- dmat[i, j]
        
        if (d < min_dist && d > 0) {
          # unit vector from j → i
          dx <- coords[i,1] - coords[j,1]
          dy <- coords[i,2] - coords[j,2]
          u <- c(dx, dy) / d
          
          # required separation move
          push <- (min_dist - d) / 2 * step
          
          coords[i,] <- coords[i,] + u * push
          coords[j,] <- coords[j,] - u * push
          
          moved <- TRUE
        }
      }
    }
    
    if (!moved) break
  }
  
  # return updated sf object
  cu_lab <- sf::st_as_sf(data.frame(x = coords[,1], y = coords[,2]),
                         coords = c("x","y"),
                         crs = sf::st_crs(points_sf))$geometry
  return(cu_lab)
}

cu_lab <- enforce_min_distance(cu_centroid, min_dist = 0.35, step = 0.05)

#------------------------------------------------------------------------------
# Read in dataframe summaries of exposure output 
#------------------------------------------------------------------------------
source("docs/online-supplement-functions.R")

expsoure <- exposure_R1 # Use exposure post-review, excluding SSS

ind <- which(exposure$rcp == rcp & exposure$period == "mid")
var <- which(names(exposure) == "exp_simple")
exp_summary0 <- cbind(
  tapply(exposure[ind, var], exposure$cuid[ind], median),
  tapply(exposure[ind, var], exposure$cuid[ind], min),
  tapply(exposure[ind, var], exposure$cuid[ind], max)
)

# Do all CUs have shapefiles?
as.numeric(row.names(exp_summary0))[which(as.numeric(row.names(exp_summary0)) %in% cu_centroid$CUID == FALSE)]

exp_summary <- exp_summary0[match(cu_order, as.numeric(row.names(exp_summary0))), ]
cu_list.plot <- cu_list[match(cu_order, cu_list$cuid),]


#------------------------------------------------------------------------------
# Colour palettes
#------------------------------------------------------------------------------
# NO species colours
quartz(width = 1.5, height = 3, pointsize = 12)
par(mar = c(0,0,0,0))
plot(rep(0.4, 6), c(6:1), pch = species_pch, bg = "white", col = 1, lwd = 2, xlim = c(0, 2), ylim = c(0, 7), bty = "n", xaxt = "n", yaxt = "n", cex = 3, xlab = "", ylab = "")
text(0.8, c(6:1), names(species_cols), adj = 0)

# Create palette for overall exposure, over range of median for rcp and period
col_levels_overall <- seq(min(exp_summary[,1]), max(exp_summary[,1]), length.out = 100)

quartz(width = 5, height = 2, pointsize = 12)
par(mar = c(3,1,4,1))
plot(0,0, "n", xlim = c(min(exp_summary[,1]), max(exp_summary[,1])), ylim = c(0, 1), bty = "n", yaxt = "n", xlab = "", ylab = "")
segments(x0 = col_levels_overall, x1 = col_levels_overall, y0 = 0, y1 = 1, col = col_palette, lwd = 5, lend = 2)

#------------------------------------------------------------------------------
# Map
#------------------------------------------------------------------------------
png(filename = "output/figures/map/map_fraser_R1.png", width = 9.5, height = 5, units = "in", res = 300, pointsize = 8)
par(mar = rep(0, 4))
# Set up blank plot
plot(st_geometry(shoreline), col = "white",  bg = "#D8EAEE", border = grey(0.6), xlim = c(-136, -114.26), ylim = c(48.56, 55.6))


# Add Fraser basin
plot(st_geometry(cu_union), border = NA, col = "#00000020", add = TRUE)

# coarser lakes and rivers
plot(st_geometry(lakes0), border = grey(0.6), lwd = 0.6, col = "#D8EAEE", add = TRUE)
plot(st_geometry(rivers0), col = grey(0.6), lwd = 0.6, add = TRUE)

# Pink points!
plot(st_geometry(cu_lab[which(cu_order == 710),]), pch = 23, col = NA, bg = col_palette[findInterval(exp_summary[which(cu_order == 710), 1], col_levels_overall)], cex = 3, add = TRUE, xpd = NA)

plot(st_geometry(cu_lab), pch = species_pch[cu_list.plot$species_pooled], bg = col_palette[findInterval(exp_summary[, 1], col_levels_overall)], cex = 3, lwd = 2, add = TRUE, xpd = NA) # col = species_cols[cu_list.plot$species_pooled], 


# Add numbers
centroid_coords <- st_coordinates(cu_lab)
white_cuid <- as.numeric(names(which(exp_summary[, 1] < 0.3 | exp_summary[, 1] > 0.45)))

text(x = centroid_coords[, "X"], y = centroid_coords[, "Y"], labels = match(cu_order, cu_centroid$CUID), cex = 0.7, font = 2, col = ifelse(cu_centroid$CUID %in% white_cuid, "#FFFFFF", "#000000"))

dev.off()
#-------
# Scale bar
# Calculate approximate scale bar length (example for a projected CRS)
# This requires knowing your map's extent and desired scale bar length in map units.
# For unprojected (lat/lon) data, this is more complex and less accurate.
u <- par('usr')
x_min <- u[1]
y_min <- u[3]
x_max <- u[2]
y_max <- u[4]

# # Define start and end points for the scale bar
# scale_bar_start_x <- u[2] - (u[2] - x_min) * 0.2 # 10% from left
# scale_bar_end_x <- scale_bar_start_x + 3 # Example: 200 km length
# scale_bar_y <- y_min + (u[4] - y_min) * 0.05 # 5% from bottom
# 
# # Draw the scale bar segments
# segments(scale_bar_start_x, scale_bar_y, scale_bar_end_x, scale_bar_y, lwd = 2)
# segments(scale_bar_start_x, scale_bar_y - 0.1, scale_bar_start_x, scale_bar_y + 0.1, lwd = 2)
# segments(scale_bar_end_x, scale_bar_y - 0.1, scale_bar_end_x, scale_bar_y + 0.1, lwd = 2)
# 
# # Add text label for the scale bar
# text(mean(c(scale_bar_start_x, scale_bar_end_x)), scale_bar_y + 0.01, "200 km", pos = 3)


#---------- Inset

png(filename = "output/figures/map/map_inset.png", width =4, height = 4, units = "in", res = 300, pointsize = 8)
par(mar = rep(0, 4))

plot(st_geometry(shoreline), col = "white",  bg = paste0(cols[5], 30), border = grey(0.6), xlim = c(-140, -116), ylim = c(42, 65))

# plot(st_geometry(BC), col = NA, border = 1, add = TRUE)
# Add Fraser basin
plot(st_geometry(cu_union), border = NA, col = "#00000020", add = TRUE)

# # Draw the scale bar segments
# segments(x_min, y_min, x_max, y_min, lwd = 2)
# segments(x_min, y_min, x_min, y_max, lwd = 2)
# segments(x_max, y_min, x_max, y_max, lwd = 2)
# segments(x_min, y_max, x_max, y_max, lwd = 2)
dev.off()