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

source("code/functions.R")

###############################################################################
# Define space and time variables for Fraser basin
###############################################################################

# Create time variable
# Note: PCIC output runs from 19450101 to 20991231, but is subsetted to date >= 1970-01-01 which is the earliest date used in the CCVAs (start of historical period)
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
grid_points <- read.csv("output/PCIC-grid-points_Fraser.csv")

# Grid cells corresponding to each CU and stage for Fraser SEL
incl.stages <- readRDS("output/freshwater-grid_included-cells.rds")

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

cu_list <- read.csv("data/cu_list.csv") %>% 
  subset(Species == "Lake sockeye" & Region == "Fraser" & Notes != "Extinct")

n.CUs <- length(unique(cu_list$Conservation.Unit))

#------------------------------------------------------------------------------
# Timing: Spawning
#------------------------------------------------------------------------------

# ** Need to consider years in freshwater rearing
timing0 <- read.csv("data/timing/4Life_cycle_table_by_CU.csv") %>% subset(species == "SEL" & region == "fraser")
timing0$cuid <- cu_list$CUID[match(timing0$culabel, cu_list$Conservation.Unit)]
timing0$culabel[is.na(timing0$cuid)] # CUs missing CUID are extinct
timing0 <- timing0[!is.na(timing0$cuid), ]


# Create new dataframe to store timing information
timing <- data.frame(
  culabel = rep(timing0$culabel, each = 4),
  cuid = rep(timing0$cuid, each = 4),
  stage = rep(c("eggs_alevin", "fw_rearing", "adult_migration", "spawning"), n.CUs),
  start = NA,
  n.days = NA # Use number of days instead of end, to account for stages that last >1 year
)

# Add timing start/end data
timing$start[timing$stage == "eggs_alevin"] <- round(timing0$inc_start)
timing$n.days[timing$stage == "eggs_alevin"] <- ifelse(
  timing0$inc_end < timing0$inc_start, # If end is earlier in the year than start (most species)
  (365 - round(timing0$inc_start) + 1) + round(timing0$inc_end),
  round(timing0$mig_end) - round(timing0$mig_start) + 1) # else Steelhead incubate over summer

# ** Need to consider years in freshwater rearing
# (All Fraser steelhead are 2+)
timing0$rear_days <- 0
timing0$rear_days[which(timing0$oe_age == "1+")] <- 365
timing0$rear_days[which(timing0$oe_age == "2+")] <- 2*365

# Ocean entry DOY is NOT ALWAYS greater than fry migration DOY 
# timing0$oe_end > timing0$fm_start
timing$start[timing$stage == "fw_rearing"] <- round(timing0$fm_start)
timing$n.days[timing$stage == "fw_rearing"] <- timing0$rear_days + round(timing0$oe_end) - round(timing0$fm_start) 

# Need to account for steelhead, who's migration may extend over winter...
timing$start[timing$stage == "adult_migration"] <- round(timing0$mig_start)
timing$n.days[timing$stage == "adult_migration"] <- ifelse(
  timing0$mig_end < timing0$mig_start, # If end is earlier in the year than start (Steelhead)
  (365 - round(timing0$mig_start) + 1) + round(timing0$mig_end),
  round(timing0$mig_end) - round(timing0$mig_start) + 1)
  
# Spawning - need to account for over winter
timing$start[timing$stage == "spawning"] <- round(timing0$sp_start)
timing$n.days[timing$stage == "spawning"] <- ifelse(
  timing0$sp_end < timing0$sp_start, # If end is earlier in the year than start (Steelhead)
  (365 - round(timing0$sp_start) + 1) + round(timing0$sp_end),
  round(timing0$sp_end) - round(timing0$sp_start) + 1)
# write.csv(timing, "output/freshwater_timing_FraserSEL.csv", row.names = FALSE)

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
                   dim = c(n.models, n.CUs, length(stages), 4, 4, 3),
                   dimnames = list(
                     gcms$modelName,
                     cu_list$Conservation.Unit, 
                     stages,
                     c("optimalTemp", "criticalTemp", "optimalFlow", "criticalFlow"),
                     c("hist", "early", "mid", "late"), 
                     c("median", "lower", "upper")))



# Create list to store full spatial output for plotting maps
# Later will calculate median across all GCMs
fw_spat <- list(); length(fw_spat) <- n.CUs*length(stages)
dim(fw_spat) <- c( n.CUs, length(stages))
dimnames(fw_spat) <- list(cu_list$Conservation.Unit, stages)

for(i in 1:n.CUs){
  for(j in 1:length(stages)){
    fw_spat[[i, j]] <- array(data = NA,
                             dim = c(n.models, 4, 4, length(incl.stages[[i,j]])),
                             dimnames = list(
                               gcms$modelName,
                               c("optimalTemp", "criticalTemp", "optimalFlow", "criticalFlow"),
                               c("hist", "early", "mid", "late"),
                               incl.stages[[i,j]])
    )
  }
}


# Store temperature thresholds for each stage
tempThresh <- data.frame(
  cuid = rep(cu_list$CUID, each = length(stages)*n.models),
  culabel = rep(cu_list$Conservation.Unit, each = length(stages)*n.models),
  stage = rep(rep(stages, each = n.models), n.CUs),
  model = rep(gcms$modelName, length(stages)*n.CUs),
  critical_max = rep(NA, n.models*n.CUs*length(stages)),
  optimal_max = rep(NA, n.models*n.CUs*length(stages)))

# Track time to run (takes ~ 20 mins for Fraser SEL CUs)
start.time <- Sys.time()
lap.time <- Sys.time()
time.for.models <- matrix(NA, nrow = n.CUs, ncol = n.models)

#------------------------------------------------------------------------------
# Loop through all Global Climate Models
#------------------------------------------------------------------------------

for(m in 1:n.models){
  
  # Load model output; this takes a minute
  GCM_var <- readRDS(paste0("data/processed-data/PCIC_", gcms$modelName[m], "_processed.rds")) 
  
  Ts.weeklyMax <- GCM_var[[1]]
  Qs.weeklyMean <- GCM_var[[2]]
  
  # Check dimnames: this is how output will be organized
  # range(as.numeric(dimnames(Ts.weeklyMax)[[1]])) # space, grid_points$id
  # range(as.Date(dimnames(Ts.weeklyMax)[[2]])) # time, date
  # Extract grid_points$id for rows of Ts.weeklyMax and Qs.weeklyMean
  # These numbers are what are referred to in incl.stages
  grid.ref <- as.numeric(dimnames(Ts.weeklyMax)[[1]])

  rm(GCM_var)
  
  # Plots
  if(3 == 2){
    
    
  }
  
  #------------------------------------------------------------------------------
  # Loop through each CU
  #------------------------------------------------------------------------------
  
  for(i in 1:n.CUs){
    
    cuid <- cu_list$CUID[i]
    
    #------------------------------------------------------------------------------
    # Loop through each life stage
    #------------------------------------------------------------------------------
    
    for(j in 1:length(stages)){
      
      # If there are timing data
      if(length(which(timing$cuid == cuid & timing$stage == stages[j])) == 0){
        warning(paste0("No timing data for ", cu_list$Conservation.Unit[i], " ", stages[j]))
        
      } else {
        
        # Extract timing and covert to day-of-year
        timing.ij <- timing[which(timing$cuid == cuid & timing$stage == stages[j])[1], c("start", "n.days")]
        
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
        
        for(p in 1:4){
          for(y in 1:30){
            start.ind <- which(period == c("hist", "early", "mid", "late")[p] & years == period.years[y, p] & DOY == timing.ij$start)
            
            # FOr the last year in late century, can't extend if start + n.days > 365
            if((start.ind + timing.ij$n.days - 1) > dim(Ts.weeklyMax)[2]){
              
              n.available <- length(start.ind:(dim(Ts.weeklyMax)[2]))
              
              Ts.ij[, p, 1:n.available, y] <- Ts.weeklyMax[match(incl, grid.ref), start.ind:(dim(Ts.weeklyMax)[2])]
              Qs.ij[, p, 1:n.available, y] <- Qs.weeklyMean[match(incl, grid.ref), start.ind:(dim(Ts.weeklyMax)[2])]
              
            } else {
              
              Ts.ij[, p, , y] <- Ts.weeklyMax[match(incl, grid.ref), start.ind:(start.ind + timing.ij$n.days - 1)]
              Qs.ij[, p, , y] <- Qs.weeklyMean[match(incl, grid.ref), start.ind:(start.ind + timing.ij$n.days - 1)]
            }
            
        }}
        
        # Plot flow
        if(3 == 2){
          
          
          
        }
        
        # Plot temperature
        if(3 == 2){
          xx <- 80 # which grid cell
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
        
        # Calculate acute (critical) and chronic (optimal) temperature thresholds
        # from historical period
        
        thresh <- apply(Ts.ij[, 1, ], 1, quantile, c(0.1, 0.9, 0.025, 0.975), na.rm = TRUE)
        rownames(thresh) <- c("opt_min", "opt_max", "crit_min", "crit_max")
        
        if(3 == 2){
          par(mfrow = c(1,1), mar = c(4,4,2,1), oma = c(0,0,2,0))
          hist(histTij[1, ], las = 1, xlab = "Weekly max. (˚C)", main = "", col = paste0(cols[5], 60), border = NA)
          abline(v = thresh[c("opt_min", "opt_max"), 10], col = cols[2], lty = 2)
          abline(v = thresh[c("crit_min", "crit_max"), 10], col = cols[2])
          for(k in 2:4){
            
            hist(Ts.weeklyMax[k, 10, ], las = 1, xlab = "", main = "", col = paste0(cols[c(3,4,1)[k-1]], 60), border = NA, add = TRUE)
          }
          mtext(side = 3, outer = TRUE, paste(cu_list$culabel[i], stages[j]))
          legend("topright", fill = paste0(cols[c(5,3,4,1)], 60), border = NA, legend = c("hist", "early", "mid", "late"))
        } # end plot
        
        #----------------------------------------------------------------------
        # Apply hard thresholds where historical temps are extreme
        #----------------------------------------------------------------------
        # ** Might need to think about also imposing higher crit_max for northern CUs where historical temp distribution is cold...
        thresh["crit_max", which(thresh["crit_max", ] > 26)] <- 26 # Critical max 26˚C, Richter & Kolmes (2005)
        thresh["opt_max", which(thresh["opt_max", ] > 24)] <- 24 # Chronic max 24˚C, Sullivan 2000
        
        #  ** Ignore minimum thresholds for now **
        # thresh["opt_min", which(thresh["opt_min", ] < c(8, 10, 8, 7)[j])] <- c(8, 10, 8, 7)[j] # Example of minimum temps from Mayer (in review)
        thresh["opt_min", ] <- 0
        thresh["crit_min", ] <- 0
        
        # Check that optimal is always inside critical
        if(length(which(thresh["crit_max", ] < thresh["opt_max", ])) > 0){
          stop("crit_max < opt_max")
        }
        # if(length(which(thresh["crit_min", ] > thresh["opt_min", ])) > 0){
        #   stop("crit_min > opt_min")
        # }
        
        #----------------------------------------------------------------------
        # Store median thresholds for output
        #----------------------------------------------------------------------
        tempThresh[which(tempThresh$cuid == cuid & tempThresh$stage == stages[j] & tempThresh$model == gcms$modelName[m]), c("critical_min", "critical_max")] <- apply(thresh[c("crit_min", "crit_max"), ], 1, median)
        
        tempThresh[which(tempThresh$cuid == cuid & tempThresh$stage == stages[j] & tempThresh$model == gcms$modelName[m]), c("optimal_min", "optimal_max")] <- apply(thresh[c("opt_min", "opt_max"), ], 1, median)
        
        #----------------------------------------------------------------------
        # Calculate days above threshold (DAT) and % DAT
        #----------------------------------------------------------------------
        Ts.dat <- Ts.perc <- array(NA, 
                        dim = c(2, 4, n.grid),
                        dimnames = list(c("opt", "crit"), c("hist", "early", "mid", 'late'), NULL)
        )
        
        for(p in 1:4){ # For each period
          for(k in 1:2){ # For opt and crit
            
            Ts.dat[k, p, ] <- apply(Ts.ij[, p, ] > matrix(rep(thresh[c(2, 4)[k], ], 30*timing.ij$n.days), n.grid, 30*timing.ij$n.days), 1, sum, na.rm = TRUE)
          
            # Percentage of days;
            # Need to account for late-century period not being able to extend into next year for 2099 if timing.ij$start + timing.ij$n.days > 365
            Ts.perc[k, p, ] <- Ts.dat[k, p, ] / sum(!is.na(Ts.ij[1, p, ]))
            } # end k
          } # end p
        
        # Store full spatial output
        fw_spat[[i,j]][m, c("optimalTemp", "criticalTemp"), , ] <- Ts.perc
        
        # store summary output
        for(p in 1:4){ # for each period
          for(k in 1:2){ # for optimal and critical
            # Store output in end array
            fw_output[m, i, j, c("optimalTemp", "criticalTemp")[k], p, ] <- c(
              median(Ts.perc[k, p, ]),
              range(Ts.perc[k, p, ]))
            
          }}
        
        #-----------------------------------------------------------------------------
        # Flow
        #-----------------------------------------------------------------------------
        
        # Calculate MAD over historical period for each included grid cell
        # (Use all DOY, not Qs.ij for life stage)
        mad <- apply(Qs.weeklyMean[match(incl, grid.ref), which(period == "hist")], 1, mean, na.rm = TRUE)
        
        # Calculate acute (critical) and chronic (optimal) temperature thresholds
        thresh_mad <- mad %*% matrix(c(0.2, 0.1), nrow = 1)
        colnames(thresh_mad) <- c("opt", "crit")
        
        # Calculate days below % mad threshold
        Qs.dbt <- Qs.perc <- array(NA, 
                        dim = c(2, 4, n.grid),
                        dimnames = list(c("opt", "crit"), 
                                        c("hist", "early", "mid", 'late'), 
                                        NULL)
        )
        for(p in 1:4){
          for(k in 1:2){
            
            Qs.dbt[k, p, ] <- apply(Qs.ij[, p, ] < matrix(rep(thresh_mad[, k], 30*timing.ij$n.days), n.grid, 30*timing.ij$n.days), 1, sum, na.rm = TRUE)
            
            # Percentage of days;
            # Need to account for late-century period not being able to extend into next year for 2099 if timing.ij$start + timing.ij$n.days > 365
            Qs.perc[k, p, ] <- Qs.dbt[k, p, ] / sum(!is.na(Qs.ij[1, p, ]))
            
            
            # Qs.dbt[k, p, ] <- calcDBT(x = Qs.weeklyMean[p, , ], 
            #                           mad_threshold = c(thresh_mad[, c("opt", "crit")[k]]))
          }}
        
        
       # Store full spatial output
        fw_spat[[i,j]][m, c("optimalFlow", "criticalFlow"), , ] <- Qs.perc
        
        for(p in 1:4){ # for each period
          for(k in 1:2){
            # Store output in end array
            fw_output[m, i, j, c("optimalFlow", "criticalFlow")[k], p, ] <- c(
              median(Qs.perc[c("opt", "crit")[k], p, ]),
              range(Qs.perc[c("opt", "crit")[k], p, ]))
          }}
        
      } # end if timing data
      # print(".")
      
    } # end life stage j
    
    print(paste0("Done ", cu_list$Conservation.Unit[i]))
    
  } # end CU i
  
  print(paste0("Done ", gcms$modelName[m]))
  
} # end model m

print(Sys.time() - start.time)
# Time difference of 7 minutes
###############################################################################
# Save output
###############################################################################

write.csv(tempThresh, file = paste0("output/tempThresholds_FraserSEL_", Sys.Date(), ".csv"))
saveRDS(fw_spat, file = paste0("output/freshwater_spat_FraserSEL_", Sys.Date(), ".rds"))
saveRDS(fw_output, file = paste0("output/freshwater_output_FraserSEL_",  Sys.Date(), ".rds"))

# write.csv(grid_points, file = "output/PCIC_grid-points_Fraser.csv")
