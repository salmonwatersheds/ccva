###############################################################################
# Code to summarize gridded hydrologic model output from the Pacific Climate 
# Impacts Consortium and calculate projected changes in exposure to stream
# temperature and flow outside of the optimal ranges for the given species 
# and life stage.
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

source("freshwater/code/freshwater-functions.R")

###############################################################################
# Define space and time variables for Fraser basin
###############################################################################

# Create time variable
# Note: PCIC output runs from 19450101 to 20991231, but is subsetted in 
# 1b_process-PCIC_fraser.R to date >= 1970-01-01 which is the earliest date used 
# in the CCVAs (start of historical period)
date <- as.Date(0:47481, origin = "1970-01-01")

DOY <- as.numeric(strftime(date, format = "%j"))
months <- as.numeric(strftime(date, format = "%m"))
years <- as.numeric(strftime(date, format = "%Y"))

# Define period variable
period <- rep(NA, length(date))
period[which(years %in% c(1970:1999))] <- "hist"
period[which(years %in% c(2010:2039))] <- "early"
period[which(years %in% c(2040:2069))] <- "mid"
period[which(years %in% c(2070:2099))] <- "late"

period.years <- cbind(
  hist = c(1970:1999), 
  early = c(2010:2039), 
  mid = c(2040:2069), 
  late = c(2070:2099))

# Spatial grid
grid_points <- read.csv("freshwater/data/processed-data/PCIC-grid-points_fraser.csv")

# Grid cells corresponding to each CU and stage for Fraser SEL
incl.stages <- readRDS("freshwater/output/PCIC_incl.rds")
incl.stages_na.rm <- incl.stages

#-----------------------------------------------------------------------------
# Define freshwater stages
#-----------------------------------------------------------------------------

stages <- c("adult_migration", "spawning", "eggs_alevin", "fw_rearing")

###############################################################################
# Population sensitivity data
###############################################################################

#------------------------------------------------------------------------------
# List of CUs in the PSE
#------------------------------------------------------------------------------

cu_list <- read.csv("data/conservationunits_decoder.csv") %>% 
  subset(region == "Fraser" & (cu_type != "Extinct" | is.na(cu_type)))

# Check that CUs equal what's in incl.stages
cu_list$pooledcuid %in% as.numeric(dimnames(incl.stages)[[1]])
as.numeric(dimnames(incl.stages)[[1]]) %in% cu_list$pooledcuid
dim(cu_list)[1] == dim(incl.stages)[1]

# Create vector of cuid
cuid <- as.numeric(dimnames(incl.stages)[[1]])
n.CUs <- length(cuid)

# Order CU list to match incl.stages order (CK, CM, CO, PKE, PKO, SEL, SER, SH)
cu_list <- cu_list[match(cuid, cu_list$cuid), ]

# Create pooled species field (for matching to temp max)
cu_list$species_pooled <- cu_list$species_name
cu_list$species_pooled[cu_list$species_name %in% c("Lake sockeye", "River sockeye")] <- "Sockeye"
cu_list$species_pooled[cu_list$species_name %in% c("Pink (odd)", "Pink (even)")] <- "Pink"

#------------------------------------------------------------------------------
# Timing
#------------------------------------------------------------------------------
# Previous versions of timing data are in data/timing/archive
# The current version has shortened adult migration duration
# for SEL CUs, using time it takes to migrate from ocean to lake and assuming
# that sockeye can then seek thermal refuge at depth in lakes.

timing <- read.csv("data/timing/timing-fraser.csv")

#------------------------------------------------------------------------------
# Temperature thresholds
#------------------------------------------------------------------------------
# If using absolute thresholds, read in from Provincial guidelines
topt <- read.csv("freshwater/data/optimum-temps_BC2001.csv")

# Adjust migration max (and add for steelhead, which is missing from BC Guidelines)
topt$max[topt$stage == "adult_migration"] <- 18

# Create matrix for max by CU
tmax <- matrix(data = NA, nrow = n.CUs, ncol = 4, dimnames = list(cuid, stages))
for(i in 1:n.CUs){
  tmax[i, ] <- topt$max[which(topt$species == cu_list$species_pooled[i])][match(stages, unique(topt$stage))]
}

###############################################################################
# Global climate model names
###############################################################################
# Vector of global climate model names (to extract netcdf output)
gcms <- data.frame(
  modelName = c(
    "ACCESS1-0",
    "CanESM2",
    "CCSM4",
    "CNRM-CM5",
    "HadGEM2-ES",
    "MPI-ESM-LR"),
  simulation = c(
    "r1i1p1",
    "r1i1p1",
    "r2i1p1",
    "r1i1p1",
    "r1i1p1",
    "r3i1p1"
  )
)

n.models <- length(gcms$modelName)

###############################################################################
# Loop through emissions scenarios, models, life stages, CUs
###############################################################################

for(r in 1:2){
  
  rcp <- c(45, 85)[r]
  
  #------------------------------------------------------------------------------
  # Setup
  #------------------------------------------------------------------------------
  
  fw_output <- array(NA,
                     dim = c(n.models, n.CUs, length(stages), 2, 4, 3),
                     dimnames = list(
                       gcms$modelName,
                       cuid, 
                       stages,
                       c("temp", "flow"),
                       c("hist", "early", "mid", "late"), 
                       c("median", "lower", "upper")))
  
  
  
  # Create list to store full spatial output for plotting maps
  # Later will calculate median across all GCMs
  fw_spat <- list(); length(fw_spat) <- n.CUs*length(stages)
  dim(fw_spat) <- c(n.CUs, length(stages))
  dimnames(fw_spat) <- list(cuid, stages)
  
  for(i in 1:n.CUs){
    for(j in 1:length(stages)){
      fw_spat[[i, j]] <- array(data = NA,
                               dim = c(n.models, 2, 4, length(incl.stages[[i,j]])),
                               dimnames = list(
                                 gcms$modelName,
                                 c("temp", "flow"),
                                 c("hist", "early", "mid", "late"),
                                 incl.stages[[i,j]])
      )
    }
  }
  
  
  # Track time to run (takes ~ 20 mins for Fraser SEL CUs)
  start.time <- Sys.time()
  lap.time <- Sys.time()
  time.for.models <- matrix(NA, nrow = n.CUs, ncol = n.models)
  
  #------------------------------------------------------------------------------
  # Loop through all Global Climate Models
  #------------------------------------------------------------------------------
  
  for(m in 1:n.models){
    
    # Load model output; this takes a minute
    GCM_var <- readRDS(paste0("freshwater/data/processed-data/PCIC_", gcms$modelName[m], "_rcp", rcp,  "_processed.rds")) 
    
    Ts.weeklyMax <- GCM_var[[1]]
    # Note that max(NA, na.rm = TRUE) return -Inf; need to set these to NA for weeklyMax output
    Ts.weeklyMax[which(Ts.weeklyMax == -Inf)] <- NA
    
    Qs.weeklyMean <- GCM_var[[2]]
    
    # Check dimnames: this is how output will be organized
    # range(as.numeric(dimnames(Ts.weeklyMax)[[1]])) # space, grid_points$id
    # range(as.Date(dimnames(Ts.weeklyMax)[[2]])) # time, date
    # Extract grid_points$id for rows of Ts.weeklyMax and Qs.weeklyMean
    # These numbers are what are referred to in incl.stages
    grid.ref <- as.numeric(dimnames(Ts.weeklyMax)[[1]])
    
    rm(GCM_var)
    
    #------------------------------------------------------------------------------
    # Loop through each CU
    #------------------------------------------------------------------------------
    
    for(i in 1:n.CUs){
      # i <- 47 #Bowron
      # i <- grep("Bay Winter", cu_list$cu_name_pse)
      #------------------------------------------------------------------------------
      # Loop through each life stage
      #------------------------------------------------------------------------------
      
      for(j in 1:length(stages)){

        if(cu_list$species_pooled[i] %in% c("Pink", "Chum") & stages[j] == "fw_rearing"){
            # No freshwater rearing for pink and chum = zero exposure during this stage
            fw_output[m, i, j, , , ] <- 0
            fw_spat[[i,j]][m, , , ] <- 0
            
          } else {
            # Extract timing and covert to day-of-year
            timing.ij <- timing[which(timing$cuid == cuid[i] & timing$stage == stages[j])[1], c("start", "duration", "duration_total")]
            
            #-----------------------------------------------------------------------------
            # Identify grid cells for given life stage
            #-----------------------------------------------------------------------------
            
            # See 1a_spatial-distribution_fraser.R for code that identifies which grid cells should be used
            # for each life stage.
            incl0 <- incl.stages[[i,j]]
            
            # In some cases, there are projections NA for all years for a grid cell (e.g., cells on the border)
            # Identify and remove these?
            if(length(incl0) == 1){ # if there's just one grid cell
              if(sum(is.na(Ts.weeklyMax[match(incl0, grid.ref),])) == ncol(Ts.weeklyMax)){
                stop("No data for single grid cell relevant for this CU and stage!")
              } else {
                incl <- incl0
              }
            } else { # if there is more than one grd cell
              no.data <- which(apply(is.na(Ts.weeklyMax[match(incl0, grid.ref),]), 1, sum) == ncol(Ts.weeklyMax))
              if(length(no.data) > 0){
                # Check that flow data is also missing
                no.flow <- which(apply(is.na(Qs.weeklyMean[match(incl0, grid.ref),]), 1, sum) == ncol(Qs.weeklyMean))
                if(length(no.data) == length(no.flow) & sum(no.data - no.flow) == 0){
                  incl <- incl0[-no.data]
                } else {
                  stop("Missing temp data not the same as missing flow data for incl.")
                }
              } else {
                incl <- incl0
              }
            }
            
            # Track which grid cells are actually used.
            incl.stages_na.rm[[i,j]] <- incl
            
            # Dimensions
            n.grid <- length(incl)
            
            #-----------------------------------------------------------------------------
            # Subset stream temp relevant to CU and stage
            #-----------------------------------------------------------------------------
            Ts.ij <- Qs.ij <- array(NA,
                                    dim = c(n.grid, 4, timing.ij$duration, 30),
                                    dimnames = list(incl, c("hist", "early", "mid", "late"), NULL, NULL))
            # Note: for stages that span >365 days, some temps may be counted twice, because essentially two cohorts occupy that space at the same time.
            
            for(p in 1:4){ # for each period
              for(y in 1:30){ # for each year in the period
                start.ind <- which(period == c("hist", "early", "mid", "late")[p] & years == period.years[y, p] & DOY == timing.ij$start)
                
                # For the last year in late century, can't extend if start + n.days > 365
                if((start.ind + timing.ij$duration - 1) > dim(Ts.weeklyMax)[2]){
                  
                  n.available <- length(start.ind:(dim(Ts.weeklyMax)[2]))
                  
                  Ts.ij[, p, 1:n.available, y] <- Ts.weeklyMax[match(incl, grid.ref), start.ind:(dim(Ts.weeklyMax)[2])]
                  Qs.ij[, p, 1:n.available, y] <- Qs.weeklyMean[match(incl, grid.ref), start.ind:(dim(Ts.weeklyMax)[2])]
                  
                } else {
                  
                  Ts.ij[, p, , y] <- Ts.weeklyMax[match(incl, grid.ref), start.ind:(start.ind + timing.ij$duration - 1)]
                  Qs.ij[, p, , y] <- Qs.weeklyMean[match(incl, grid.ref), start.ind:(start.ind + timing.ij$duration - 1)]
                }
                
              } # end year y
            } # end period p
            
            # PLOTTING
            if(3 == 2){
              # Plot flow
              xDate <- as.Date(paste("1998", 1:365, sep = "-"), format = "%Y-%j")
              xDate_display <- as.Date(c(paste("1998", 1:12, rep(1, 12), sep = "-"), "1999-01-01"), format = "%Y-%m-%d")
              
              ht <- which(period == c("hist"))
              er <- which(period == c("early"))
              md <- which(period == c("mid"))
              lt <- which(period == c("late"))
              xx <- 10
              # xx <- 198
              quartz(width = 8, height = 5, pointsize = 12)
              plot(xDate, tapply(Qs.weeklyMean[match(incl, grid.ref)[xx], ht], DOY[ht], median)[1:365], "l", col = cols[5], lwd = 2, log = "y", xaxt = "n", xaxs = "i", xlab = "", ylab = "Flow (cm/s)")
              axis(side = 1, at = xDate_display, labels = c(month.abb, month.abb[1]))
              abline(v = xDate_display, lty = 3, col = grey(0.8), lwd = 0.8)
              
              abline(h = 0.2*mean(Qs.weeklyMean[match(incl, grid.ref)[xx], ht], na.rm = TRUE), lty = 2)
              
              polygon(x = c(xDate, rev(xDate)),
                      y = c(tapply(Qs.weeklyMean[match(incl, grid.ref)[xx], ht], DOY[ht], quantile, 0.25)[1:365],
                            rev(tapply(Qs.weeklyMean[match(incl, grid.ref)[xx], ht], DOY[ht], quantile, 0.75)[1:365])),
                      col = paste0(cols[5], 30),
                      border = NA)
              
              lines(xDate, tapply(Qs.weeklyMean[match(incl, grid.ref)[xx], md], DOY[md], median)[1:365], lwd = 2, col = cols[3])
              polygon(x = c(xDate, rev(xDate)),
                      y = c(tapply(Qs.weeklyMean[match(incl, grid.ref)[xx], md], DOY[md], quantile, 0.25)[1:365],
                            rev(tapply(Qs.weeklyMean[match(incl, grid.ref)[xx], md], DOY[md], quantile, 0.75)[1:365])),
                      col = paste0(cols[3], 30),
                      border = NA)
              
              # # Late century
              # lines(xDate, tapply(Qs.weeklyMean[match(incl, grid.ref)[xx], l], DOY[l], median)[1:365], lwd = 2, col = cols[1])
              # polygon(x = c(xDate, rev(xDate)),
              #         y = c(tapply(Qs.weeklyMean[match(incl, grid.ref)[xx], l], DOY[l], quantile, 0.25)[1:365],
              #               rev(tapply(Qs.weeklyMean[match(incl, grid.ref)[xx], l], DOY[l], quantile, 0.75)[1:365])),
              #         col = paste0(cols[1], 30),
              #         border = NA)
              # --------
              # Stream temp
              plot(xDate, tapply(Ts.weeklyMax[match(incl, grid.ref)[xx], h], DOY[ht], median)[1:365], "l", col = cols[5], lwd = 2, ylim = c(0, 25))
              
              polygon(x = c(xDate, rev(xDate)),
                      y = c(tapply(Ts.weeklyMax[match(incl, grid.ref)[xx], h], DOY[ht], quantile, 0.25)[1:365],
                            rev(tapply(Ts.weeklyMax[match(incl, grid.ref)[xx], h], DOY[ht], quantile, 0.75)[1:365])),
                      col = paste0(cols[5], 30),
                      border = NA)
              
              lines(xDate, tapply(Ts.weeklyMax[match(incl, grid.ref)[xx], m], DOY[md], median)[1:365], lwd = 2, col = cols[3])
              
              polygon(x = c(xDate, rev(xDate)),
                      y = c(tapply(Ts.weeklyMax[match(incl, grid.ref)[xx], m], DOY[md], quantile, 0.25)[1:365],
                            rev(tapply(Ts.weeklyMax[match(incl, grid.ref)[xx], m], DOY[md], quantile, 0.75)[1:365])),
                      col = paste0(cols[3], 30),
                      border = NA)
              
              # Plot propotion of days
              quartz(width = 12, height = 3, pointsize = 12)
              par(mar = c(4,4,2,1))
              plot(date[seq(1, 47482, 14)], Ts.weeklyMax[match(incl, grid.ref)[xx], seq(1, 47482, 14)], "l", las = 1, ylab = "Stream temperature", xlab = "Date", bty = "l")
              polygon(x = date[c(range(h), rev(range(h)))], y = c(-1, -1, 40, 40), border = NA, col = paste0(cols[5], 20))
              polygon(x = date[c(range(e), rev(range(e)))], y = c(-1, -1, 40, 40), border = NA, col = paste0(cols[3], 20))
              polygon(x = date[c(range(m), rev(range(m)))], y = c(-1, -1, 40, 40), border = NA, col = paste0(cols[4], 20))
              polygon(x = date[c(range(l), rev(range(l)))], y = c(-1, -1, 40, 40), border = NA, col = paste0(cols[1], 20))
              
              # Plot days over threshold
              plot(1970:1999, apply(Ts.ij[xx, 1, , ] > tmax[i,j], 2, sum), xlim = c(1950, 2100), col = paste0(cols[5], 40), ylim = c(0, timing.ij$duration), ylab = "Number of days over threshold", xlab = "", pch = 19, bty = "l", las = 1)
              segments(x0 = 1970, x1 = 1999, y0 = sum(Ts.ij[xx, 1, , ] > tmax[i,j])/30, y1 = sum(Ts.ij[xx, 1, , ] > tmax[i,j])/30, col = cols[5], lwd = 3)
              
              for(p in 2:4){
                points(period.years[, p], apply(Ts.ij[xx, p, , ] > tmax[i,j], 2, sum), col = paste0(cols[c(3,4,1)[p-1]], 40), pch = 19)
                segments(x0 = period.years[1, p], x1 = period.years[30, p], y0 = sum(Ts.ij[xx, p, , ] > tmax[i,j])/30, y1 = sum(Ts.ij[xx, p, , ] > tmax[i,j])/30, col = cols[c(3,4,1)[p-1]], lwd = 3)
              }
              abline(v = c(1969.5, 1999.5, 2009.5, 2039.5, 2069.5), col = grey(0.8))
              
              # Proportion of stage over threshold
              plot(1970:1999, apply(Ts.ij[xx, 1, , ] > tmax[i,j], 2, sum)/timing.ij$duration, xlim = c(1950, 2100), col = paste0(cols[5], 40), ylim = c(0, 1), ylab = "Proportion of time over threshold", xlab = "", pch = 19, bty = "l", las = 1)
              segments(x0 = 1970, x1 = 1999, y0 = sum(Ts.ij[xx, 1, , ] > tmax[i,j])/30/timing.ij$duration, y1 = sum(Ts.ij[xx, 1, , ] > tmax[i,j])/30/timing.ij$duration, col = cols[5], lwd = 3)
              
              for(p in 2:4){
                points(period.years[, p], apply(Ts.ij[xx, p, , ] > tmax[i,j], 2, sum)/timing.ij$duration, col = paste0(cols[c(3,4,1)[p-1]], 40), pch = 19)
                segments(x0 = period.years[1, p], x1 = period.years[30, p], y0 = sum(Ts.ij[xx, p, , ] > tmax[i,j])/30/timing.ij$duration, y1 = sum(Ts.ij[xx, p, , ] > tmax[i,j])/30/timing.ij$duration, col = cols[c(3,4,1)[p-1]], lwd = 3)
              }
              abline(v = c(1969.5, 1999.5, 2009.5, 2039.5, 2069.5), col = grey(0.8))
              
            }
            
            # Plot temperature dots over threshold
            if(3 == 2){
              xx <- 1 # which grid cell
              start.ind <- which(period == c("hist", "early", "mid", "late")[p] & years == period.years[y, p] & DOY == timing.ij$start)
              
              # Plot 5 years
              par(mar = c(4,4,4,1), bg = "white")
              plot(date[years %in% period.years[1:5, p]], Ts.weeklyMax[match(incl[xx], grid.ref), which(years %in% period.years[1:5, p])], "l", xlab = "Year", ylab = "Weekly max. temp (˚C)", las = 1, bty = "l", xaxs ="i", col = grey(0.5))
              
              for(y in 1:5){
                start.ind <- which(years %in% period.years[y, p] & DOY == timing.ij$start)
                points(date[start.ind:(start.ind + timing.ij$duration - 1)], Ts.ij[xx, p, , y], pch = 19, cex =0.6)
                segments(x0 = date[start.ind],
                         y0 = 31 + (y-1),
                         x1 = date[(start.ind + timing.ij$duration - 1)],
                         y1 = 31 + (y-1),
                         col = paste0(cols[y], 60), 
                         lwd = 10,
                         xpd = NA
                )
              }
              mtext(side = 3, paste(cu_list$Conservation.Unit[i], cu_list$Species[i], "-", stages[j]), line = 3)
              abline(h = pmax(quantile(Ts.ij[xx, 1, , ], c(0.9, 0.975), na.rm = TRUE), c(24, 26)), lty = c(2, 3)) 
              
            } # end if 3 == 2 plot
            
            # Collapse Ts.ij years and days for easy
            dim(Ts.ij) <- c(n.grid, 4, 30*timing.ij$duration) 
            dim(Qs.ij) <- c(n.grid, 4, 30*timing.ij$duration) 
            
            #-----------------------------------------------------------------------------
            # Stream temperature 
            #-----------------------------------------------------------------------------
            
            #----------------------------------------------------------------------
            # Calculate days above threshold (DAT) and % DAT
            #----------------------------------------------------------------------
            Ts.dat <- Ts.perc <- array(NA, 
                                       dim = c(4, n.grid),
                                       dimnames = list(c("hist", "early", "mid", 'late'), NULL)
            )
            
            for(p in 1:4){ # For each period
              Ts.dat[p, ] <- apply(Ts.ij[, p, ] > matrix(rep(tmax[i, j], 30*timing.ij$duration), n.grid, 30*timing.ij$duration), 1, sum, na.rm = TRUE)
              
              # Percentage of days;
              # Need to account for late-century period not being able to extend into next year for 2099 if timing.ij$start + timing.ij$duration > 365
              # If there are missing data (Ts.ij[1, p, ] == NA), exclude those in the denominator
              # Note that missing data will be differnt for each grid cell
                  Ts.perc[p, ] <- Ts.dat[p, ] / (timing.ij$duration_total*30 - ifelse(n.grid == 1, sum(is.na(Ts.ij[, p, ])), apply(is.na(Ts.ij[, p, ]), 1, sum)))
            
              } # end p
            
            # For grid cells that didn't have any data, set to NA
            # These are grid cells at the mouth of the Fraser or in the USA
            if(length(incl) == 1){
              tmean <- mean(Ts.weeklyMax[match(incl, grid.ref), which(period == "hist")], na.rm = TRUE)
            } else {
              tmean <- apply(Ts.weeklyMax[match(incl, grid.ref), which(period == "hist")], 1, mean, na.rm = TRUE)
              }
            Ts.perc[, which(tmean == -Inf)] <- NA
            
            # Store full spatial output
            fw_spat[[i,j]][m, "temp", , match(incl, incl0)] <- Ts.perc
            
            # store summary output
            for(p in 1:4){ # for each period
              # Store output in end array
              fw_output[m, i, j, "temp", p, ] <- c(
                median(Ts.perc[p, ], na.rm = TRUE),
                range(Ts.perc[p, ], na.rm = TRUE))
              
            } # end p
            
            #-----------------------------------------------------------------------------
            # Flow
            #-----------------------------------------------------------------------------
            
            # Calculate MAD over historical period for each included grid cell
            # (Use all DOY, not Qs.ij for life stage)
            if(n.grid > 1){
              mad <- apply(Qs.weeklyMean[match(incl, grid.ref), which(period == "hist")], 1, mean, na.rm = TRUE)
            } else { # if just one grid cell
              mad <- mean(Qs.weeklyMean[match(incl, grid.ref), which(period == "hist")], na.rm = TRUE)
            }
            
            # Calculate acute (critical) and chronic (optimal) temperature thresholds
            thresh_mad <- mad %*% matrix(c(0.2, 0.1), nrow = 1)
            colnames(thresh_mad) <- c("opt", "crit")
            
            # Calculate days below % mad threshold
            Qs.dbt <- Qs.perc <- array(NA, 
                                       dim = c(4, n.grid),
                                       dimnames = list(c("hist", "early", "mid", 'late'), 
                                                       NULL)
            )
            
            k <- 1 # Use 20% MAD threshold for now
            
            for(p in 1:4){
              Qs.dbt[p, ] <- apply(Qs.ij[, p, ] < matrix(rep(thresh_mad[, k], 30*timing.ij$duration), n.grid, 30*timing.ij$duration), 1, sum, na.rm = TRUE)
              
              # Percentage of days;
              # Need to account for late-century period not being able to extend into next year for 2099 if timing.ij$start + timing.ij$duration > 365
              # Qs.perc[p, ] <- Qs.dbt[p, ] / sum(!is.na(Qs.ij[which(!is.na(mad))[1], p, ]))
              # Qs.perc[p, ] <- Qs.dbt[p, ] / ((timing.ij$duration_total - sum(is.na(Ts.ij[1, p, ]))) * 30)
              
              Qs.perc[p, ] <- Qs.dbt[p, ] / (timing.ij$duration_total*30 - ifelse(n.grid == 1, sum(is.na(Ts.ij[, p, ])), apply(is.na(Qs.ij[, p, ]), 1, sum)))
              
              
              
            } # end p
            
            # For grid cells that didn't have any data, set to NA
            # These are grid cells at the mouth of the Fraser or in the USA
            Qs.perc[, is.na(mad)] <- NA
            
            # Store full spatial output
            fw_spat[[i,j]][m, "flow", , match(incl, incl0)] <- Qs.perc
            
            for(p in 1:4){ # for each period
              # Store output in end array
              fw_output[m, i, j, "flow", p, ] <- c(
                median(Qs.perc[p, ], na.rm = TRUE),
                range(Qs.perc[p, ], na.rm = TRUE))
            }
            
          } # if if not pink/chum fw_rearing
        
      } # end life stage j
      
      rm(incl, incl0)
      print(paste0("Done ", cuid[i]))
      
    } # end CU i
    
    print(paste0("Done ", gcms$modelName[m]))
    
  } # end model m
  
  print(Sys.time() - start.time)
  # Time difference of 7 minutes
  ###############################################################################
  # Save output
  ###############################################################################
  
  saveRDS(fw_spat, file = paste0("freshwater/output/freshwater_spat_fraser_rcp", rcp, "_", Sys.Date(), "_shortSELmig.rds"))
  saveRDS(fw_output, file = paste0("freshwater/output/freshwater_output_fraser_rcp", rcp, "_",  Sys.Date(), "_shortSELmig.rds"))
  
}

saveRDS(incl.stages_na.rm, file = "freshwater/output/PCIC_incl_narm.rds")
saveRDS(incl.stages_na.rm, file = paste0("freshwater/output/PCIC_incl_narm_",  Sys.Date(), ".rds"))
