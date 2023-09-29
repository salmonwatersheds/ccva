###############################################################################
# Code to summarize projected changes in exposure to stream temperature and 
# flow outside of the historical norms for the given species and life stage.
# Takes output from calc-freshwater-exposure.R and tries to make sense of it!
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

setwd("freshwater/")
###############################################################################
# Load background spatial layers
###############################################################################

# Set broad bounds for clipping broad spatial layers
bounds <- c(xmin = -124.5, ymin = 50.85, xmax = -123.6, ymax = 51.7)

#------------------------------------------------------------------------------
# Lakes, rivers, and shorelines
#------------------------------------------------------------------------------

lakes <- readRDS("data/spatial/layers/waterbodies_250.rds")# %>% st_crop(bounds)
rivers <- readRDS("data/spatial/layers/watercourse_250.rds") #%>% st_crop(bounds)

bounds0 <- c(xmin = -127, ymin = 49, xmax = -116.5, ymax = 56)
lakes0 <- readRDS("data/spatial/layers/waterbodies_lowRes.rds")# %>% st_crop(bounds)
rivers0 <- readRDS("data/spatial/layers/watercourse_lowRes.rds") #%>% st_crop(bounds)
BC <- readRDS("data/spatial/layers/BC_lowRes.rds")
shoreline <- st_read(dsn = "data/spatial/layers/GSHHS_i_L1.shp")

#------------------------------------------------------------------------------
# Spawning Zone of Influence: Fraser SEL
#------------------------------------------------------------------------------

spawn_zoi <- st_read(dsn = "data/spatial/ZOI/fraser-spawning-zoi/fraser_spawning_zoi_SEL.shp") %>% st_transform(crs = 4269)

#------------------------------------------------------------------------------
# Adult migration lines: Fraser SEL
#------------------------------------------------------------------------------
mig_paths <- st_read(dsn = "data/spatial/fw-migration/spawn_timing_migration_paths_wCU_fraser.shp") %>% st_transform(crs = 4269)

#------------------------------------------------------------------------------
# Conservation Unit Boundaries: SEL
#------------------------------------------------------------------------------

# DFO shapefiels from Open Data Canada
# shp_path <- "data/spatial/CU-boundaries/Lake_Type_Sockeye_Salmon_CU_Boundary/SEL_CU_BOUNDARY_En"
# cu_boundary <- st_read(dsn = paste0(shp_path, ".shp"), promote_to_multi = FALSE)

# PSF version of shapefiles
shp_path <- "data/spatial/CU-boundaries/PSF_CUs_updatedDec2021"
cu_boundary <- st_read(dsn = paste0(shp_path, ".shp"), promote_to_multi = FALSE) %>%
  subset(regionname == "Fraser" & species %in% c("Lake sockeye", "Sockeye-Lake"))

###############################################################################
# Create space variables for plotting
###############################################################################

grid_points <- read.csv("output/PCIC-grid-points_Fraser.csv") %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4269)

#------------------------------------------------------------------------------
# List of CUs in the PSE
#------------------------------------------------------------------------------

cu_list <- read.csv("../data/cu_list.csv") %>% 
  subset(Species == "Lake sockeye" & Region == "Fraser" & Notes != "Extinct")

n.CUs <- length(unique(cu_list$Conservation.Unit))

#------------------------------------------------------------------------------
# Timing
#------------------------------------------------------------------------------

timing <- read.csv("output/freshwater_timing_FraserSEL.csv")
###############################################################################
# CU Overlay plots:
# Calculate median and range across models for each CU, by life-stage
###############################################################################

fw_output <- readRDS("output/freshwater_output_FraserSEL_2023-09-28.rds")
# dimensions:
#   gcms$modelName,
#   cu_list$Conservation.Unit, 
#   stages,
#   c("optimalTemp", "criticalTemp", "optimalFlow", "criticalFlow"),
#   c("hist", "early", "mid", "late"), 
#   c("median", "lower", "upper")))


# Define stages
stages <- dimnames(fw_output)[[3]]
stages.all <- c(stages, "early_marine", "marine_rearing")
stage.names <- c(c("Adult\nmigration", "Spawning", "Incubation", "Freshwater\nrearing", "Early\nmarine", "Marine\nrearing"))
oj <- rev(c(3, 4, 5, 6, 1, 2)) # Order in which we want stages to appear

# How much time does each CU spend in each stage?
numDays <- array(NA, dim = c(n.CUs, 4), dimnames = list(cu_list$Conservation.Unit, stages))
for(i in 1:n.CUs){
  for(j in 1:4){
    numDays[i,j] <- timing$n.days[which(timing$culabel == cu_list$Conservation.Unit[i] & timing$stage == stages[j])]
  }}

# # Create dataframe to store CU summary outputs
# fw_cu <- data.frame(
#   region = rep("Fraser", n.CUs * length(stages) * 4),
#   species = rep("Lake sockeye", n.CUs * length(stages) * 4),
#   Conservation.Unit = rep(cu_list$Conservation.Unit, each = length(stages) * 4),
#   cuid = rep(cu_list$CUID, each = length(stages) * 4),
#   stage = rep(rep(stages, each = 4), n.CUs),
#   period = rep(dimnames(fw_output)[[5]])
# )

# set colour palette
col_levels <- seq(0, 0.9, 0.1)
n <- length(col_levels)
col_palette <- paste0(colorRampPalette(colors = cols[c(5,3,1)])(n = length(col_levels) + 1))

pdf(file = paste0("output/figures/FW-exposure_CUlevel_FraserSEL_", Sys.time(), ".pdf"), width = 6.5, height = 6)
# Loop through each CU and plot; or select CU and variable to plot
for(i in 1:n.CUs){
  for(k in 1:4){
    
# i <- 1
# k <- 1
    
    par(mar =c(4, 8, 5, 2), bg = "white")
    plot(c(0, 1), c(0.5, 6.5), "n", yaxt ="n", bty = "n", xlab = "Proportion of days", ylab = "", xaxs = "i", yaxs = "i", xaxt = "n")
    
    # Add colour legend at bottom
    for(q in 1:n){
      polygon(x = c(0, 1/n, 1/n, 0)+ (q-1) * 1/n, y = c(0, 0, 0.5, 0.5), border = NA, col = col_palette[q], xpd = NA)
    }
    
    axis(side = 1, tck = -0.02)
    abline(h = c(1:6) + 0.5, xpd = NA)
    segments(x0 = 0, y0 = 0.5, x1 = 0, y1 = 6.5, xpd = NA)
    segments(x0 = seq(0, 1, 0.1), y0 = 0.5, x1 = seq(0, 1, 0.1), y1 = 6.5, col = grey(0.8), xpd = NA)
    segments(x0 = 0, y0 = 0.5, x1 = -10, y1 = 0.5, xpd = NA)
    
    text(rep(-0.1, 6), c(1:6), stage.names[oj], xpd = NA)
    
    # Grey out marine stages
    polygon(x = c(-10, 10, 10, -10), y = c(2.5, 2.5, 4.5, 4.5), col = "#00000030", border = NA, xpd = NA)
    text(rep(0.5, 2), c(3, 4), "Not Applicable", col = "white", font = 2)
    
    for(j in 1:4){
      for(p in 1:4){
        x <- median(fw_output[, i, j, k, p, 1]) 
        pCol <- col_palette[findInterval(x, col_levels)]
        
        points(
          x, 
          oj[j] + c(-0.3, -0.1, 0.1, 0.3)[p], 
          col = pCol,
          bg = pCol,
          pch = c(22, 25, 21, 24)[p],
          cex = 1.5, xpd = NA)
        segments(
          x0 = min(fw_output[, i, j, k, p, 1]),
          y0 = oj[j] + c(-0.3, -0.1, 0.1, 0.3)[p],
          x1 = max(fw_output[, i, j, k, p, 1]),
          y1 = oj[j] + c(-0.3, -0.1, 0.1, 0.3)[p],
          col = pCol,
          lwd = 1.5
        )
        
      } # end p
      text(1.03, oj[j], numDays[i,j], xpd = NA, font = 2)
    } # end j
    
    legend(-0.2, 7.2, ncol = 4, pch = c(22, 25, 21, 24), pt.cex = 1.5, col = 1, pt.bg = 1, c("Historical", "Early-century", "Mid-century", "Late-century"), xpd = NA, bty = "n", bg = NA)
    
    # mtext(side = 3, adj = 1, c("Optimal  Critical"), line = 2)
    # mtext(side = 3, adj = 1, c("Optimal","Critical")[k], line = 3.2)
    
    mtext(side = 3, paste0("Fraser SEL - ", cu_list$Conservation.Unit[i], "\n", c("Optimal Stream Temperature","Critical Stream Temperature", "Optimal Flow", "Critical Flow")[k]), line = 3)
  
    } # end k
} # end cus
dev.off()

###############################################################################
# CU map:
# Show spatial grid of exposure for a CU for mid-century period only, 
# with life-stages overlaid 
###############################################################################

# Set colour palette for fw exposure (0 - 1 % of days)
col_levels <- seq(0, 0.9, 0.1)
col_palette <- paste0(colorRampPalette(colors = cols[c(5,3,1)])(n = length(col_levels) + 1), 80)
col_palette_FULL <- colorRampPalette(colors = cols[c(5,3,1)])(n = length(col_levels) + 1)

# Load grid cell ids for each CU and life stage
incl.stages <- readRDS("output/freshwater-grid_included-cells.rds")

# Use mid-century period (can't show all periods overlaid on the map...)
p <- 3

# Select CU (or loop through)

i <- 19
cuid <- cu_list$CUID[i]

# Subset spatial data for selected CU

if(length(which(cu_boundary$FULL_CU_IN == cu_list$Full.CU.Index[i])) == 0){
  stop(paste0("No CU boundary for ", cu_list$Conservation.Unit[i]))
} else {
  cu_boundary.i <- cu_boundary[which(cu_boundary$FULL_CU_IN == cu_list$Full.CU.Index[i]),]
}

zoi.i <- spawn_zoi[which(spawn_zoi$cuid == cuid), ]

mig_paths.i <- mig_paths[which(mig_paths$cuid == cuid), ]

# Set bounds for map as extent of migration and spawning ZOI
bounds0 <- st_bbox(mig_paths.i)[c(1,3,2,4)]
bounds1 <- st_bbox(zoi.i)[c(1,3,2,4)]
bounds2 <- apply(cbind(bounds0, bounds1), 1, range)
bounds <- c(min(bounds2[,1]), max(bounds2[,2]), min(bounds2[,3]), max(bounds2[,4]))

# Plot map
par(bg = "white")

# Set up blank plot
plot(st_geometry(BC), border = NA, col = "white", axes = FALSE, las = 1, ylim = bounds[3:4], xlim = bounds[1:2], bty = "o")

# Add shoreline
plot(shoreline, add = TRUE, col = NA, border = 1)

# coarser lakes and rivers
plot(st_geometry(lakes0), border = 1, col = NA, add = TRUE)
plot(st_geometry(rivers0), col = 1, add = TRUE)

# # Add CU ZOI and boundary?
# if(!is.na(zoi.i$cuid)){
#   plot(st_geometry(zoi.i), border = grey(0.6), col = grey(0.9), add = TRUE)
# }
# 
# if(!is.na(cu_boundary.i$CU_NAME)){
#   plot(st_geometry(cu_boundary.i), border = grey(0.4), col = grey(0.7), lwd = 1, add = TRUE)
# }

d <- 1/16

# Plot grid cells, coloured by exposure
for(j in 1:4){
  
  # Highlight those cells
  incl <- incl.stages[[i,j]]
  grid_polys_incl <- st_as_sf(data.frame(
    id = rep(grid_points$id[incl], each = 5),
    rep(c("SW0", "NW", "NE", "SE", "SW1"), length(incl)),
    lon = c(rep(rep(lon, length(lat))[incl], each = 5) + rep(c(- d/2, -d/2, d/2,  d/2, -d/2), length(incl))),
    lat = c(rep(rep(lat, each = length(lon))[incl], each = 5) + rep(c(- d/2,  d/2, d/2, -d/2, -d/2), length(incl)))
  ), coords = c("lon", "lat"), crs = 4269) %>% 
    group_by(id) %>%
    summarise(geometry = st_combine(geometry)) %>%
    st_cast("POLYGON") 


exp_level <- apply(fw_spat[[i,j]][, k, p, ], 2, median, na.rm = TRUE) # For some reason, spatial output only came through for m = 6...

col_polys <- col_palette[findInterval(exp_level, col_levels)]
plot(grid_polys_incl, col = col_polys, border = col_polys, add = TRUE)

}

mtext(side = 3, cu_list$Conservation.Unit[i])
#---# 
# Highlight a lifestage
#---#
j <- 3
d <- 1/16
incl <- incl.stages[[i,j]]
grid_polys_incl <- st_as_sf(data.frame(
  id = rep(grid_points$id[incl], each = 5),
  rep(c("SW0", "NW", "NE", "SE", "SW1"), length(incl)),
  lon = c(rep(rep(lon, length(lat))[incl], each = 5) + rep(c(- d/2, -d/2, d/2,  d/2, -d/2), length(incl))),
  lat = c(rep(rep(lat, each = length(lon))[incl], each = 5) + rep(c(- d/2,  d/2, d/2, -d/2, -d/2), length(incl)))
), coords = c("lon", "lat"), crs = 4269) %>% 
  group_by(id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") 

exp_level <- apply(fw_spat[[i,j]][, k, p, ], 2, median, na.rm = TRUE)
col_polys <- col_palette_FULL[findInterval(exp_level, col_levels)]

plot(grid_polys_incl, col = col_polys, border = col_polys, add = TRUE)
plot(st_union(grid_polys_incl), col = NA, border = 1, lwd = 1.5, add = TRUE)
plot(grid_polys_incl[20,], col = 1, border = 1, add = TRUE)

mtext(side = 3, line = 2, stages[j], font = 2)


###############################################################################
# Full species view: compare CUs
###############################################################################

# Take the weighted average of pdays across stages, median across models
fw_summary <- array(NA, dim = c(n.CUs, 4, 4, 3), dimnames = list(cu_list$Conservation.Unit, c("optimalTemp",  "criticalTemp", "optimalFlow",  "criticalFlow"), c("hist",  "early", "mid",   "late"), c("median", "lower",  "upper" )))
for(i in 1:n.CUs){
  for(k in 1:4){
    for(p in 1:4){
  
      fw_summary[i, k, p, 1] <- sum(apply(fw_output[, i, , k, p, 1], 2, median) * numDays[i, ])
      fw_summary[i, k, p, 2] <- sum(apply(fw_output[, i, , k, p, 1], 2, min) * numDays[i, ])
      fw_summary[i, k, p, 3] <- sum(apply(fw_output[, i, , k, p, 1], 2, max) * numDays[i, ])
}}}

col_levels <- seq(min(fw_summary), max(fw_summary), length.out = 10)
n <- length(col_levels)
col_palette <- paste0(colorRampPalette(colors = cols[c(5,3,1)])(n = length(col_levels) + 1))


pdf(file = "output/figures/FW-exposure_FullSpecies_FraserSEL.pdf", width = 8.5, height = 6)
for(k in 1:4){
  oCU <- order(fw_summary[, k, 4, 1], decreasing = FALSE)
  
  par(mar = c(4, 20, 2, 1))
  plot(range(fw_summary), c(0.5, n.CUs+0.5), "n", bty = "n", yaxt = "n", ylab = "", yaxs = "i", xaxs = "i", xaxt = "n", xlab = "Total days outside optimal")
  u <- par('usr')
  text(rep(-0.02*(u[2] - u[1]), n.CUs), c(1:n.CUs), cu_list$Conservation.Unit[oCU], xpd = NA, adj = 1, cex = 0.8)
  
  for(q in 1:n){
    polygon(x = c(0, max(fw_summary)/n, max(fw_summary)/n, 0)+ (q-1) * max(fw_summary)/n, y = c(0, 0, 0.5, 0.5), border = NA, col = col_palette[q], xpd = NA)
  }
  
  axis(side = 1, tck = -0.02)
  abline(h = c(1:n.CUs)+0.5, xpd = NA, col = grey(0.8))
  segments(x0 = 0, y0 = 0.5, x1 = 0, y1 = n.CUs+0.5, xpd = NA)
  segments(x0 = 0, y0 = 0.5, x1 = -10, y1 = 0.5, xpd = NA)
  
  for(i in 1:n.CUs){
    for(p in 1:4){
      x <- fw_summary[oCU[i], k, p, 1]
      pCol <- col_palette[findInterval(x, col_levels)]
      
      points(
        x, 
        i + c(-0.3, -0.1, 0.1, 0.3)[p], 
        col = pCol,
        bg = pCol,
        pch = c(22, 25, 21, 24)[p],
        cex = 1, xpd = NA)
      segments(
        x0 = fw_summary[oCU[i], k, p, 2],
        y0 = i + c(-0.3, -0.1, 0.1, 0.3)[p],
        x1 = fw_summary[oCU[i], k, p, 3],
        y1 = i + c(-0.3, -0.1, 0.1, 0.3)[p],
        col = pCol
      )
      
    } # end p
  } # end j
  mtext(side = 3, paste0("Fraser SEL - ", c("Optimal Temperature","Critical Temperature", "Optimal Flow", "Critical Flow")[k]))
  
} #End k

dev.off()


######### Map of CUs

bounds <- st_bbox(cu_boundary)[c(1,3,2,4)]

par(bg = "white")
plot(st_geometry(BC), border = NA, col = NA, axes = FALSE, las = 1, ylim = bounds[3:4], xlim = bounds[1:2], bty = "o")
plot(BC, add = TRUE, col = NA, border = 1)

# coarser lakes and rivers
plot(st_geometry(lakes0), border = 1, col = NA, add = TRUE)
plot(st_geometry(rivers0), col = 1, add = TRUE)

p <- 4

for(i in 1:n.CUs){
  
  x <- fw_summary[i, k, p, 1]
  pCol <- col_palette[findInterval(x, col_levels)]
  
  plot(st_geometry(cu_boundary[which(cu_boundary$CUID == cu_list$CUID[i]), ]), col = paste0(pCol, 30), border = pCol, add = TRUE)
  
}


###############################################################################
###############################################################################
# Extra Plots
###############################################################################
###############################################################################


###############################################################################
# Flow: why does prop days below threshold go down under late-century?
###############################################################################

dim(Qs.ij) # Keep in long-form; e.g., [n.gird, periods, n.days, year]

Qs.ij3 <- Qs.ij
dim(Qs.ij3) <- c(n.grid, 4, timing.ij$n.days*30)

mad <- apply(Qs.weeklyMean[match(incl, grid.ref), which(period == "hist")], 1, mean, na.rm = TRUE)
xDate <- as.Date(paste(2000, c(1:366), sep = "-"), format = "%Y-%j")

Qs <- Qs.weeklyMean[match(incl, grid.ref), ]
Qs.med <- Qs[20,]#apply(Qs, 2, median)

plot(range(xDate), c(0, 7), "n", bty = "l")

polygon(x = xDate[c(timing.ij$start, 366, 366, timing.ij$start)], y = c(-1, -1, 8, 8), col = grey(0.8), border = NA)
polygon(x = xDate[c(1, rep(timing.ij$n.days + timing.ij$start - 365, 2), 1)], y = c(-1, -1, 8, 8), col = grey(0.8), border = NA)


for(p in 1:4){
  lines(xDate, tapply(Qs.med[which(period == c("hist", "early", "mid", "late")[p])], DOY[which(period == c("hist", "early", "mid", "late")[p])], median), col = cols[c(5,3,4,1)[p]], lwd = 2, xpd = NA)
}

abline(h = 0.2*mad[20])


###############################################################################
# Summarize output for Shiny app
###############################################################################

#-----------------------------------------------------------------------------
# Spatial summary: median across models
#-----------------------------------------------------------------------------

# Median across models for each grid cell, each CU, each stage; mid-century ONLY
fw_spat_summary <- list()
length(fw_spat_summary) <- n.CUs*length(stages) 
dim(fw_spat_summary) <- c(n.CUs, length(stages))

p <- 3 # Mid century only

for(i in 1:n.CUs){
  for(j in 1:4){
    fw_spat_summary[[i,j]] <- array(data = NA, 
                                    dim = c(length(incl.stages[[i,j]]), 2))
    for(k in 1:2){
      fw_spat_summary[[i,j]][, k] <- apply(fw_spat[[i,j]][, c(1,3)[k], p, ], 2, median, na.rm = TRUE)
    }
  }
}
   
saveRDS(fw_spat_summary, "docs/data/output/fw_spat_summary.rds")                                                   
#-----------------------------------------------------------------------------
# Output summary: median and range across models
#-----------------------------------------------------------------------------
fw_output_summary <- array(NA, dim = c(n.CUs, length(stages), 2, 4, 3),
                           dimnames = list(cu_list$Conservation.Unit, stages, c("Stream temperature", "Low flow"), c("hist", "early", "mid", "late"), c("median", "min", "max")))

for(i in 1:n.CUs){
  for(j in 1:4){
    for(k in 1:2){
      for(p in 1:4){
        fw_output_summary[i, j, k, p, 1] <- median(fw_output[, i, j, c(1,3)[k], p, 1])
        fw_output_summary[i, j, k, p, 2] <- min(fw_output[, i, j, c(1,3)[k], p, 1])
        fw_output_summary[i, j, k, p, 3] <- max(fw_output[, i, j, c(1,3)[k], p, 1])
      }
    }
  }
}

saveRDS(fw_output_summary, "docs/data/output/fw_output_summary.rds")                                              