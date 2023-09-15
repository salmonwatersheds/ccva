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

###############################################################################
# Create space variables for plotting
###############################################################################

grid_points <- read.csv("output/PCIC_grid-points_Fraser.csv") %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4269)

#------------------------------------------------------------------------------
# List of CUs in the PSE
#------------------------------------------------------------------------------

cu_list <- read.csv("data/cu_list.csv") %>% 
  subset(Species == "Lake sockeye" & Region == "Fraser" & Notes != "Extinct")

n.CUs <- length(unique(cu_list$Conservation.Unit))


###############################################################################
# Calculate median and range across models for each CU and life-stage
###############################################################################

fw_output <- readRDS("output/freshwater_output_FraserSEL.rds")
# dimensions:
#   gcms$modelName,
#   cu_list$Conservation.Unit, 
#   stages,
#   c("optimalTemp", "criticalTemp", "optimalFlow", "criticalFlow"),
#   c("hist", "early", "mid", "late"), 
#   c("median", "lower", "upper")))

stages <- dimnames(fw_output)[[3]]

fw_ccva <- data.frame(
  region = rep("Fraser", n.CUs * length(stages) * 4),
  species = rep("Lake sockeye", n.CUs * length(stages) * 4),
  Conservation.Unit = rep(cu_list$Conservation.Unit, each = length(stages) * 4),
  cuid = rep(cu_list$CUID, each = length(stages) * 4),
  stage = rep(rep(stages, each = 4), n.CUs),
  period = rep(dimnames(fw_output)[[5]])
)

# Example for a given CU
i <- 4 # Chilko early summer
j <- 1

stages.all <- c(stages, "early_marine", "marine_rearing")
oj <- rev(c(3, 4, 5, 6, 1, 2))
stage.names <- c(c("Adult\nmigration", "Spawning", "Incubation", "Freshwater\nrearing", "Early\nmarine", "Marine\nrearing"))

pCols <- cols[c(5,3,4,1)]
# Stream temperature
fw_output[, i, j, 1:2, , 1]


par(mar =c(4, 5, 2, 1))
plot(c(0, 1), c(0.5, 6.5), "n", yaxt ="n", bty = "n", xlab = "Proportion of days", ylab = "", xaxs = "i", yaxs = "i")

abline(h = c(1:6)+0.5, xpd = NA)
segments(x0 = 0, y0 = 0.5, x1 = 0, y1 = 6.5, xpd = NA)
segments(x0 = 0, y0 = 0.5, x1 = -10, y1 = 0.5, xpd = NA)

text(rep(-0.1, 6), c(1:6), stage.names[oj], xpd = NA)

# Grey out marine stages
polygon(x = c(-10, 10, 10, -10), y = c(2.5, 2.5, 4.5, 4.5), col = "#00000030", border = NA, xpd = NA)
k <- 2

for(j in 1:4){
  for(p in 1:4){
     points(
        median(fw_output[, i, j, k, p, 1]), 
        oj[j] + c(-0.3, -0.1, 0.1, 0.3)[p], 
        col = pCols[p],
        pch = 19,
        cex = 1.5)
      segments(
        x0 = min(fw_output[, i, j, k, p, 1]),
        y0 = oj[j] + c(-0.3, -0.1, 0.1, 0.3)[p],
        x1 = max(fw_output[, i, j, k, p, 1]),
        y1 = oj[j] + c(-0.3, -0.1, 0.1, 0.3)[p],
        col = pCols[p],
        lwd = 1.5
      )

  } # end p
} # end j

legend(0, 7.2, ncol = 4, pch = 19, pt.cex = 1.5, col = pCols, c("Historical", "Early-century", "Mid-century", "Late-century"), xpd = NA, bty = "n", bg = NA)

mtext(side = 3, adj = 1, c("Optimal  Critical"), line = 1.5)
mtext(side = 3, adj = 1, c("Optimal","Critical")[k], line = 3)


########
# Spatial part
########

#-------------------------------------------
# Look at grid polys for each life stage
#-----------------------------------------------------------------------------

# List of total c
incl.stages <- list(); length(incl.stages) <- n.CUs*4; dim(incl.stages) <- c(n.CUs, 4)
for(i in 1:n.CUs){
  
  cuid <- cu_list$CUID[i]
  
  # Subset spatial data for selected CU
  
  if(length(which(cu_boundary$FULL_CU_IN == cu_list$Full.CU.Index[i])) == 0){
    stop(paste0("No CU boundary for ", cu_list$Conservation.Unit[i]))
  } else {
    cu_boundary.i <- cu_boundary[which(cu_boundary$FULL_CU_IN == cu_list$Full.CU.Index[i]),]
  }
  
  zoi.i <- spawn_zoi[which(spawn_zoi$cuid == cuid), ]
  
  mig_paths.i <- mig_paths[which(mig_paths$cuid == cuid), ]
  
  for(j in 1:4){
    # For adult migration and freshwater rearing, include downstream migration route
    if(stages[j] %in% c("adult_migration", "fw_rearing")){
      intrscts <- st_is_within_distance(grid_points, mig_paths.i, 3000, sparse = FALSE)  
      incl.stages[[i,j]] <- which(c(intrscts[, 1]) == TRUE)
      # There may be more than one line segment
      if(dim(intrscts)[2] > 1){
        for(k in 2:dim(intrscts)[2]){
          incl.stages[[i,j]] <- unique(c(incl.stages[[i,j]], which(c(intrscts[, k]) == TRUE)))
        }
      }
    } 
    
    # For spawning, eggs_alevin, use spawning ZOI
    if(stages[j] %in% c("spawning", "eggs_alevin")){
      intrscts <- st_intersects(grid_points, zoi.i, sparse = FALSE)
      incl.stages[[i,j]] <- which(c(intrscts) == TRUE)
    }
    
    # For freshwater rearing, include rearing lake (CU boundary)
    if(stages[j] == "fw_rearing"){
      
      # CU boundary is a multipolygon and was throwing error: 
      p <- st_make_valid(cu_boundary.i)
      
      intrscts <- st_intersects(grid_points, p, sparse = FALSE)
      incl.stages[[i,j]] <- unique(c(incl.stages[[i,j]], which(c(intrscts) == TRUE)))
    } 
    
    # Exclude cells that have no data
    incl.stages[[i,j]] <- incl.stages[[i,j]][apply(Ts[incl.stages[[i,j]], ], 1, function(x) sum(!is.na(x))) > 0]
  }
  
} # end i

# Highlight those cells
k <- 2
p <- 4

# Map
bounds <- st_bbox(mig_paths.i)[c(1,3,2,4)]
bounds2 <- st_bbox(zoi.i)[c(1,3,2,4)]

bounds3 <- apply(cbind(bounds, bounds2), 1, range)
bounds <- c(min(bounds3[,1]), max(bounds3[,2]), min(bounds3[,3]), max(bounds3[,4]))

par(bg = "white")
plot(st_geometry(BC), border = NA, col = NA, axes = FALSE, las = 1, ylim = bounds[3:4], xlim = bounds[1:2], bty = "o")
plot(BC, add = TRUE, col = NA, border = 1)

if(!is.na(zoi.i$cuid)){
  plot(st_geometry(zoi.i), border = 1, col = NA, lwd = 1.5, add = TRUE)
}

if(!is.na(cu_boundary.i$CU_NAME)){
  plot(st_geometry(cu_boundary.i), border = 1, col = NA, lwd = 1, add = TRUE)
}

# coarser lakes and rivers
plot(st_geometry(lakes0), border = 1, col = NA, add = TRUE)
plot(st_geometry(rivers0), col = 1, add = TRUE)


# Plot cells
for(j in 1:4){
grid_polys_incl <- st_as_sf(data.frame(
  id = rep(grid_points$id[incl.stages[[i,j]]], each = 5),
  rep(c("SW0", "NW", "NE", "SE", "SW1"), length(incl.stages[[i,j]])),
  lon = c(rep(rep(lon, length(lat))[incl.stages[[i,j]]], each = 5) + rep(c(- d/2, -d/2, d/2,  d/2, -d/2), length(incl.stages[[i,j]]))),
  lat = c(rep(rep(lat, each = length(lon))[incl.stages[[i,j]]], each = 5) + rep(c(- d/2,  d/2, d/2, -d/2, -d/2), length(incl.stages[[i,j]])))
), coords = c("lon", "lat"), crs = 4269) %>% 
  group_by(id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") 


col_levels <- seq(0, 0.9, 0.1)
col_palette <- paste0(colorRampPalette(colors = cols[c(5,3,1)])(n = length(col_levels) + 1), 80)

exp_level <- apply(fw_spat[[i,j]][, k, p, ], 2, median, na.rm = TRUE) # For some reason, spatial output only came through for m = 6...

col_polys <- col_palette[findInterval(exp_level, col_levels)]
plot(grid_polys_incl, col = col_polys, border = col_polys, add = TRUE)
}

plot(grid_polys_incl, col = col_polys, border = 1, lwd = 0.5, add = TRUE)


##############
# Different visualization

col_levels <- seq(0, 0.9, 0.1)
n <- length(col_levels)
col_palette <- paste0(colorRampPalette(colors = cols[c(5,3,1)])(n = length(col_levels) + 1))


pdf("output/figures/CU_StreamTemp.pdf", width = 6, height = 5.5)
pdf("output/figures/CU_LowFlow.pdf", width = 6, height = 5.5)

for(i in 1:n.CUs){
  for(k in 1:2){

par(mar =c(4, 5, 4, 1))
plot(c(0, 1), c(0.5, 6.5), "n", yaxt ="n", bty = "n", xlab = "Proportion of days", ylab = "", xaxs = "i", yaxs = "i", xaxt = "n")

for(q in 1:n){
  polygon(x = c(0, 1/n, 1/n, 0)+ (q-1) * 1/n, y = c(0, 0, 0.5, 0.5), border = NA, col = col_palette[q], xpd = NA)
}

axis(side = 1, tck = -0.02)
abline(h = c(1:6)+0.5, xpd = NA)
segments(x0 = 0, y0 = 0.5, x1 = 0, y1 = 6.5, xpd = NA)
segments(x0 = 0, y0 = 0.5, x1 = -10, y1 = 0.5, xpd = NA)

text(rep(-0.1, 6), c(1:6), stage.names[oj], xpd = NA)

# Grey out marine stages
polygon(x = c(-10, 10, 10, -10), y = c(2.5, 2.5, 4.5, 4.5), col = "#00000030", border = NA, xpd = NA)

for(j in 1:4){
  for(p in 1:4){
    x <- median(fw_output[, i, j, c(3,4)[k], p, 1]) # Low flow
    x <- median(fw_output[, i, j, k, p, 1]) # Stream temp
    pCol <- col_palette[findInterval(x, col_levels)]
    
    points(
      x, 
      oj[j] + c(-0.3, -0.1, 0.1, 0.3)[p], 
      col = pCol,
      bg = pCol,
      pch = c(22, 25, 21, 24)[p],
      cex = 1.5)
    segments(
      x0 = min(fw_output[, i, j, k, p, 1]),
      y0 = oj[j] + c(-0.3, -0.1, 0.1, 0.3)[p],
      x1 = max(fw_output[, i, j, k, p, 1]),
      y1 = oj[j] + c(-0.3, -0.1, 0.1, 0.3)[p],
      col = pCol,
      lwd = 1.5
    )
    
  } # end p
} # end j

legend(0, 7.2, ncol = 4, pch = c(22, 25, 21, 24), pt.cex = 1.5, col = 1, pt.bg = 1, c("Historical", "Early-century", "Mid-century", "Late-century"), xpd = NA, bty = "n", bg = NA)

# mtext(side = 3, adj = 1, c("Optimal  Critical"), line = 2)
# mtext(side = 3, adj = 1, c("Optimal","Critical")[k], line = 3.2)

mtext(side = 3, paste0("Fraser SEL - ", cu_list$Conservation.Unit[i], " - ", c("Optimal","Critical")[k], " Stream Temperature"), line = 2.5)
} # end k
} # end cus
dev.off()


##########################################################
# Compare CUs
# How much time does each CU spend in each stage?
numDays <- array(NA, dim = c(n.CUs, 4), dimnames = list(cu_list$Conservation.Unit, stages))
for(i in 1:n.CUs){
  cuid <- cu_list$CUID[i]
  for(j in 1:4){
    timing.ij <- timing[which(timing$cuid == cuid & timing$stage == stages[j])[1], c("start", "end")] %>% 
      as.numeric()
    if(timing.ij[2] < timing.ij[1]){
      timing.ij[2] <- timing.ij[2] + 365
    } 
  numDays[i,j] <- timing.ij[2] - timing.ij[1] + 1
  }}


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

oCU <- order(fw_summary[, 1, 4, 1], decreasing = TRUE)

par(mar = c(4, 12, 2, 1))
plot(range(fw_summary), c(0.5, n.CUs+0.5), "n", bty = "n", yaxt = "n", ylab = "", yaxs = "i", xaxs = "i", xaxt = "n", xlab = "Total days outside optimal")
text(rep(-0.02, n.CUs), c(1:n.CUs), cu_list$Conservation.Unit[oCU], xpd = NA, adj = 1)

for(q in 1:n){
  polygon(x = c(0, max(fw_summary)/n, max(fw_summary)/n, 0)+ (q-1) * max(fw_summary)/n, y = c(0, 0, 0.5, 0.5), border = NA, col = col_palette[q], xpd = NA)
}

axis(side = 1, tck = -0.02)
abline(h = c(1:n.CUs)+0.5, xpd = NA, col = grey(0.8))
segments(x0 = 0, y0 = 0.5, x1 = 0, y1 = n.CUs+0.5, xpd = NA)
segments(x0 = 0, y0 = 0.5, x1 = -10, y1 = 0.5, xpd = NA)
k <- 1
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