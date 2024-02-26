###############################################################################
# Code to summarize gridded hydrologic model output from the Pacific Climate 
# Impacts Consortium and calculate projected changes in exposure to stream
# temperature and flow outside of the historical norms for the given species 
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

# Which rcp to run? 45 or 85
# Change here manually; output labelled accordingly (instead of looping through, since
# code takes a while to run.)
rcp <- 45
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

# ** Need to consider years in freshwater rearing
timing0 <- read.csv("data/timing/4Life_cycle_table_by_CU.csv") %>% 
  subset(region == "fraser")
timing0$SQ_culabel <- paste(timing0$species, timing0$culabel, sep = "_")

# Are all CUs there?
paste(cu_list$species_abbr, cu_list$cu_name_pse, sep = "_") %in% timing0$SQ_culabel # Yes, for fraser anyway

# Subset timing data to only look at non-extinct and non-transplant CUs, match order to cu_list
# This removes 9 CUs that Sam had; 8 SEL and one CK (Fraser-Harrison Fall Transplant (Fall 4-1))
timing0 <- timing0[match(paste(cu_list$species_abbr, cu_list$cu_name_pse, sep = "_"), timing0$SQ_culabel), ]

# Create new dataframe to store timing information
timing <- data.frame(
  species = rep(timing0$species, each = 4),
  cu_name_pse = rep(timing0$culabel, each = 4),
  cuid = rep(cuid, each = 4),
  stage = rep(c("eggs_alevin", "fw_rearing", "adult_migration", "spawning"), n.CUs),
  start = NA, # DOY that each stage starts
  n.days = NA # Use number of days instead of end, to account for stages that last >1 year
)

#------------------
# Add timing start and stage n.days
#------------------
# INCUBATION
timing$start[timing$stage == "eggs_alevin"] <- round(timing0$inc_start)
timing$n.days[timing$stage == "eggs_alevin"] <- ifelse(
  timing0$inc_end < timing0$inc_start, # If end is earlier in the year than start (most species)
  (365 - round(timing0$inc_start) + 1) + round(timing0$inc_end),
  round(timing0$inc_end) - round(timing0$inc_start) + 1) # else Steelhead incubate over summer

# REARING
# ** Need to consider years in freshwater rearing
timing0$rear_days <- 0
timing0$rear_days[which(timing0$oe_age == "1+")] <- 365
timing0$rear_days[which(timing0$species == "SH")] <- 2*365 # (All Fraser steelhead are 2+)

# Pink and chum: no rearing
timing0$fm_start[which(timing0$species %in% c("CM", "PKE", "PKO"))]

timing$start[timing$stage == "fw_rearing"] <- round(timing0$fm_start)
# timing$start[timing$stage == "fw_rearing" & timing$cuid %in% cu_list$pooledcuid[cu_list$species_abbr %in% c("CM", "PKO", "PKE")]] <- round(timing0$oe_start[which(timing0$species %in% c("CM", "PKE", "PKO"))])

timing$n.days[timing$stage == "fw_rearing"] <- timing0$rear_days + round(timing0$oe_end) - timing$start[timing$stage == "fw_rearing"]

# ADULT MIGRATION
# Need to account for steelhead, who's migration may extend over winter...
timing$start[timing$stage == "adult_migration"] <- round(timing0$mig_start)
timing$n.days[timing$stage == "adult_migration"] <- ifelse(
  timing0$mig_end < timing0$mig_start, # If end is earlier in the year than start (Steelhead)
  (365 - round(timing0$mig_start) + 1) + round(timing0$mig_end),
  round(timing0$mig_end) - round(timing0$mig_start) + 1)
  
# SPAWNING
# Need to account for over winter
timing$start[timing$stage == "spawning"] <- round(timing0$sp_start)
timing$n.days[timing$stage == "spawning"] <- ifelse(
  timing0$sp_end < timing0$sp_start, # If end is earlier in the year than start (Steelhead)
  (365 - round(timing0$sp_start) + 1) + round(timing0$sp_end),
  round(timing0$sp_end) - round(timing0$sp_start) + 1)
# write.csv(timing, file = "freshwater/output/freshwater_timing_fraser.csv", row.names = FALSE)

#------------------------------------------------------------------------------
# Temperature thresholds
#------------------------------------------------------------------------------
# If using absolute thresholds, read in from Provincial guidelines
topt <- read.csv("freshwater/data/optimum-temps_BC2001.csv")

# Adjust migration max (and add for steelhead, which i smissing from BC Guidelines)
topt$max[topt$stage == "adult_migration"] <- 18

# Create matrix for max by CU
tmax <- matrix(data = NA, nrow = n.CUs, ncol = 4, dimnames = list(cuid, stages))
for(i in 1:n.CUs){
  tmax[i, ] <- topt$max[which(topt$species == cu_list$species_pooled[i])][match(stages, unique(topt$stage))]
}

# # Compare to historical thresholds for sockeye
# hist_thresh <- read.csv("freshwater/output/tempThresholds_FraserSEL_2023-09-28.csv")
# par(mfrow = c(2,2), mar = c(3,4,2,1), oma = c(2,0,2,0))
# for(j in 1:4){ # for each stage
#   hist(hist_thresh$optimal_max[hist_thresh$stage == stages[j]], breaks = seq(10.5, 24.5, 1), col = paste0(cols[1], 30), border = NA, yaxs = "i", main = stages[j], xlab = "")
#   abline(v = topt$max[topt$stage == stages[j] & topt$Species == "Sockeye"], col= cols[1], lwd = 2)
#   abline(v = 12, lty = 2, col = cols[1])
#   
#   if(j == 2){
#     legend("topright", bty = "n", col = c(paste0(cols[1], 30), cols[1], cols[1]), lwd = c(10, 2, 1), lty = c(1, 1, 2), legend = c("Historical 90th", "Prov BC guideline max", "Forced min(Historical 90th)"))
#   }
# }
# mtext(side = 1, outer = TRUE, "Temperature (˚C)")

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
# Loop through models, life stages, CUs
###############################################################################

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
    #------------------------------------------------------------------------------
    # Loop through each life stage
    #------------------------------------------------------------------------------
    
    for(j in 1:length(stages)){
      
      # If there are timing data
      if(length(which(timing$cuid == cuid[i] & timing$stage == stages[j])) == 0){
        stop(paste0("No timing data for cuid ", cu_list$cuid[i], " ", stages[j]))
        
      } else {
        
        if(cu_list$species_pooled[i] %in% c("Pink", "Chum") & stages[j] == "fw_rearing"){
          # No freshwater rearing for pink and chum = zero exposure during this stage
          fw_output[m, i, j, , , ] <- 0
          fw_spat[[i,j]][m, , , ] <- 0
          
        } else {
          # Extract timing and covert to day-of-year
          timing.ij <- timing[which(timing$cuid == cuid[i] & timing$stage == stages[j])[1], c("start", "n.days")]
          
          #-----------------------------------------------------------------------------
          # Identify grid cells for given life stage
          #-----------------------------------------------------------------------------
          
          # See freshwater-grid.R for code that identifies which grid cells should be used
          # for each life stage.
          
          incl <- incl.stages[[i,j]]
          
          # Dimensions
          n.grid <- length(incl)
          
          #-----------------------------------------------------------------------------
          # Subset stream temp relevant to CU and stage
          #-----------------------------------------------------------------------------
          Ts.ij <- Qs.ij <- array(NA,
                                  dim = c(n.grid, 4, timing.ij$n.days, 30),
                                  dimnames = list(incl, c("hist", "early", "mid", "late"), NULL, NULL))
          # Note: for stages that span >365 days, some temps may be counted twice, because essentially two cohorts occupy that space at the same time.
          
          for(p in 1:4){ # for each period
            for(y in 1:30){ # for each year in the period
              start.ind <- which(period == c("hist", "early", "mid", "late")[p] & years == period.years[y, p] & DOY == timing.ij$start)
              
              # For the last year in late century, can't extend if start + n.days > 365
              if((start.ind + timing.ij$n.days - 1) > dim(Ts.weeklyMax)[2]){
                
                n.available <- length(start.ind:(dim(Ts.weeklyMax)[2]))
                
                Ts.ij[, p, 1:n.available, y] <- Ts.weeklyMax[match(incl, grid.ref), start.ind:(dim(Ts.weeklyMax)[2])]
                Qs.ij[, p, 1:n.available, y] <- Qs.weeklyMean[match(incl, grid.ref), start.ind:(dim(Ts.weeklyMax)[2])]
                
              } else {
                
                Ts.ij[, p, , y] <- Ts.weeklyMax[match(incl, grid.ref), start.ind:(start.ind + timing.ij$n.days - 1)]
                Qs.ij[, p, , y] <- Qs.weeklyMean[match(incl, grid.ref), start.ind:(start.ind + timing.ij$n.days - 1)]
              }
              
            }}
          
          if(3 == 2){
            # Plot flow
            xDate <- as.Date(paste("1998", 1:365, sep = "-"), format = "%Y-%j")
            
            h <- which(period == c("hist"))
            e <- which(period == c("early"))
            m <- which(period == c("mid"))
            l <- which(period == c("late"))
            xx <- 10
            # xx <- 198
        quartz(width = 8, height = 5, pointsize = 12)
            plot(xDate, tapply(Qs.weeklyMean[match(incl, grid.ref)[xx], h], DOY[h], median)[1:365], "l", col = cols[5], lwd = 2, log = "y")
            abline(h = 0.2*mean(Qs.weeklyMean[match(incl, grid.ref)[xx], which(period == "hist")], na.rm = TRUE), lty = 2)
            
            polygon(x = c(xDate, rev(xDate)),
            y = c(tapply(Qs.weeklyMean[match(incl, grid.ref)[xx], h], DOY[h], quantile, 0.25)[1:365],
                  rev(tapply(Qs.weeklyMean[match(incl, grid.ref)[xx], h], DOY[h], quantile, 0.75)[1:365])),
            col = paste0(cols[5], 30),
            border = NA)
            
            lines(xDate, tapply(Qs.weeklyMean[match(incl, grid.ref)[xx], m], DOY[m], median)[1:365], lwd = 2, col = cols[3])
            
            polygon(x = c(xDate, rev(xDate)),
            y = c(tapply(Qs.weeklyMean[match(incl, grid.ref)[xx], m], DOY[m], quantile, 0.25)[1:365],
                  rev(tapply(Qs.weeklyMean[match(incl, grid.ref)[xx], m], DOY[m], quantile, 0.75)[1:365])),
            col = paste0(cols[3], 30),
            border = NA)
            
            # Stream temp
            plot(xDate, tapply(Ts.weeklyMax[match(incl, grid.ref)[xx], h], DOY[h], median)[1:365], "l", col = cols[5], lwd = 2, ylim = c(0, 25))
            
            polygon(x = c(xDate, rev(xDate)),
                    y = c(tapply(Ts.weeklyMax[match(incl, grid.ref)[xx], h], DOY[h], quantile, 0.25)[1:365],
                          rev(tapply(Ts.weeklyMax[match(incl, grid.ref)[xx], h], DOY[h], quantile, 0.75)[1:365])),
                    col = paste0(cols[5], 30),
                    border = NA)
            
            lines(xDate, tapply(Ts.weeklyMax[match(incl, grid.ref)[xx], m], DOY[m], median)[1:365], lwd = 2, col = cols[3])
            
            polygon(x = c(xDate, rev(xDate)),
                    y = c(tapply(Ts.weeklyMax[match(incl, grid.ref)[xx], m], DOY[m], quantile, 0.25)[1:365],
                          rev(tapply(Ts.weeklyMax[match(incl, grid.ref)[xx], m], DOY[m], quantile, 0.75)[1:365])),
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
              plot(1970:1999, apply(Ts.ij[xx, 1, , ] > tmax[i,j], 2, sum), xlim = c(1950, 2100), col = paste0(cols[5], 40), ylim = c(0, timing.ij$n.days), ylab = "Number of days over threshold", xlab = "", pch = 19, bty = "l", las = 1)
                segments(x0 = 1970, x1 = 1999, y0 = sum(Ts.ij[xx, 1, , ] > tmax[i,j])/30, y1 = sum(Ts.ij[xx, 1, , ] > tmax[i,j])/30, col = cols[5], lwd = 3)
 
                for(p in 2:4){
                  points(period.years[, p], apply(Ts.ij[xx, p, , ] > tmax[i,j], 2, sum), col = paste0(cols[c(3,4,1)[p-1]], 40), pch = 19)
                segments(x0 = period.years[1, p], x1 = period.years[30, p], y0 = sum(Ts.ij[xx, p, , ] > tmax[i,j])/30, y1 = sum(Ts.ij[xx, p, , ] > tmax[i,j])/30, col = cols[c(3,4,1)[p-1]], lwd = 3)
                }
                abline(v = c(1969.5, 1999.5, 2009.5, 2039.5, 2069.5), col = grey(0.8))
                
                # Proportion of stage over threshold
                plot(1970:1999, apply(Ts.ij[xx, 1, , ] > tmax[i,j], 2, sum)/timing.ij$n.days, xlim = c(1950, 2100), col = paste0(cols[5], 40), ylim = c(0, 1), ylab = "Proportion of time over threshold", xlab = "", pch = 19, bty = "l", las = 1)
                segments(x0 = 1970, x1 = 1999, y0 = sum(Ts.ij[xx, 1, , ] > tmax[i,j])/30/timing.ij$n.days, y1 = sum(Ts.ij[xx, 1, , ] > tmax[i,j])/30/timing.ij$n.days, col = cols[5], lwd = 3)
                
                for(p in 2:4){
                  points(period.years[, p], apply(Ts.ij[xx, p, , ] > tmax[i,j], 2, sum)/timing.ij$n.days, col = paste0(cols[c(3,4,1)[p-1]], 40), pch = 19)
                  segments(x0 = period.years[1, p], x1 = period.years[30, p], y0 = sum(Ts.ij[xx, p, , ] > tmax[i,j])/30/timing.ij$n.days, y1 = sum(Ts.ij[xx, p, , ] > tmax[i,j])/30/timing.ij$n.days, col = cols[c(3,4,1)[p-1]], lwd = 3)
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
              points(date[start.ind:(start.ind + timing.ij$n.days - 1)], Ts.ij[xx, p, , y], pch = 19, cex =0.6)
              segments(x0 = date[start.ind],
                       y0 = 31 + (y-1),
                       x1 = date[(start.ind + timing.ij$n.days - 1)],
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
          dim(Ts.ij) <- c(n.grid, 4, 30*timing.ij$n.days) 
          dim(Qs.ij) <- c(n.grid, 4, 30*timing.ij$n.days) 
          
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
            Ts.dat[p, ] <- apply(Ts.ij[, p, ] > matrix(rep(tmax[i, j], 30*timing.ij$n.days), n.grid, 30*timing.ij$n.days), 1, sum, na.rm = TRUE)
            
            # Percentage of days;
            # Need to account for late-century period not being able to extend into next year for 2099 if timing.ij$start + timing.ij$n.days > 365
            Ts.perc[p, ] <- Ts.dat[p, ] / sum(!is.na(Ts.ij[1, p, ]))
          } # end p
          
          
          # Store full spatial output
          fw_spat[[i,j]][m, "temp", , ] <- Ts.perc
          
          # store summary output
          for(p in 1:4){ # for each period
            # Store output in end array
            fw_output[m, i, j, "temp", p, ] <- c(
              median(Ts.perc[p, ]),
              range(Ts.perc[p, ]))
            
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
            Qs.dbt[p, ] <- apply(Qs.ij[, p, ] < matrix(rep(thresh_mad[, k], 30*timing.ij$n.days), n.grid, 30*timing.ij$n.days), 1, sum, na.rm = TRUE)
            
            # Percentage of days;
            # Need to account for late-century period not being able to extend into next year for 2099 if timing.ij$start + timing.ij$n.days > 365
            Qs.perc[p, ] <- Qs.dbt[p, ] / sum(!is.na(Qs.ij[1, p, ]))
            
          } # end p
          
          
          # Store full spatial output
          fw_spat[[i,j]][m, "flow", , ] <- Qs.perc
          
          for(p in 1:4){ # for each period
            # Store output in end array
            fw_output[m, i, j, "flow", p, ] <- c(
              median(Qs.perc[p, ]),
              range(Qs.perc[p, ]))
          }
          
        } # if if not pink/chum fw_rearing
      } # end if there are timing data
      
    } # end life stage j
    
    print(paste0("Done ", cuid[i]))
    
  } # end CU i
  
  print(paste0("Done ", gcms$modelName[m]))
  
} # end model m

print(Sys.time() - start.time)
# Time difference of 7 minutes
###############################################################################
# Save output
###############################################################################

saveRDS(fw_spat, file = paste0("freshwater/output/freshwater_spat_fraser_rcp", rcp, "_", Sys.Date(), ".rds"))
saveRDS(fw_output, file = paste0("freshwater/output/freshwater_output_fraser_rcp", rcp, "_",  Sys.Date(), ".rds"))

# write.csv(grid_points, file = "output/PCIC_grid-points_Fraser.csv")
