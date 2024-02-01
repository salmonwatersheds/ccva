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

cols <- wes_palette("Darjeeling1")

source("marine/code/marine-functions.R")

# Which rcp to run? 45 or 85
# Change here manually; output labelled accordingly (instead of looping through, since
# code takes a while to run.)
rcp <- 45

###############################################################################
# Define space and time variables 
###############################################################################

# Create time variable
date <- as.Date(paste(01, rep(1:12, 131), rep(1970:2100, each = 12), sep = "-"), format = "%d-%m-%Y")

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
spat_EM <- read.csv("marine/output/incl_NOAA_EM.csv")[, 2]
spat_MR <- read.csv("marine/output/incl_NOAA_marRear.csv")[, 2]
spat_overall <- sort(unique(c(spat_EM, spat_MR)))

#-----------------------------------------------------------------------------
# Define marine stages
#-----------------------------------------------------------------------------

stages <- c("early_marine", "marine_rearing")

###############################################################################
# Population sensitivity data
###############################################################################

#------------------------------------------------------------------------------
# List of CUs in the PSE
#------------------------------------------------------------------------------

cu_list <- read.csv("data/conservationunits_decoder.csv") %>% 
  subset(region == "Fraser" & (cu_type != "Extinct" | is.na(cu_type)))

# Create vector of cuid
cuid <- cu_list$cuid[order(cu_list$species_pooled)] # Keep same order as in freshwater
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

# Note: reorg of timing data happens in 2_calc-freshwater-exposure.R
timing <- read.csv("freshwater/output/freshwater_timing_fraser.csv")

#------------------------------------------------------------------------------
# Temperature thresholds
#------------------------------------------------------------------------------
# If using absolute thresholds, read in from Provincial guidelines
topt <- read.csv("data/optimum-temps.csv") %>%
  filter(stage == "marine")

# Create matrix for max by CU
tmax <- topt$max[match(cu_list$species_pooled, topt$species)]

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
    "MPI-ESM-LR"))

n.models <- length(gcms$modelName)

###############################################################################
# Loop through models, life stages, CUs
###############################################################################

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------

mar_output <- array(NA,
                   dim = c(n.models, n.CUs, length(stages), 2, 4, 3),
                   dimnames = list(
                     gcms$modelName,
                     cuid, 
                     stages,
                     c("SST", "SSS"),
                     c("hist", "early", "mid", "late"), 
                     c("median", "lower", "upper")))



# Create list to store full spatial output for plotting maps
# Later will calculate median across all GCMs
mar_spat <- list(); length(mar_spat) <- n.CUs*length(stages)
dim(mar_spat) <- c(n.CUs, length(stages))
dimnames(mar_spat) <- list(cuid, stages)

for(i in 1:n.CUs){
  for(j in 1:length(stages)){
    mar_spat[[i, j]] <- array(data = NA,
                             dim = c(n.models, 2, 4, length(list(spat_EM, spat_MR)[[j]])),
                             dimnames = list(
                               gcms$modelName,
                               c("SST", "SSS"),
                               c("hist", "early", "mid", "late"),
                               list(spat_EM, spat_MR)[[j]])
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
  SST <- loadNOAA(
    variable = "SST", # Which variable to load? One of "PP", "pH", "SSS", "SST"
    model = "Can-ESM2", # Which GCM? One of "ACCESS1-0", "Can-ESM2", "CCSM4", "CNRM-CM5", "HADGEM2-ES", "MPI-ESM-LR"
    rcp = 45, # Which emissions scenario? Options 45 or 85
    spat_id = spat_overall 
  )
  
  SSS <- loadNOAA(
    variable = "SST", # Which variable to load? One of "PP", "pH", "SSS", "SST"
    model = "Can-ESM2", # Which GCM? One of "ACCESS1-0", "Can-ESM2", "CCSM4", "CNRM-CM5", "HADGEM2-ES", "MPI-ESM-LR"
    rcp = 45, # Which emissions scenario? Options 45 or 85
    spat_id = spat_overall 
  )
  
  
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
    
    #------------------------------------------------------------------------------
    # Loop through each life stage
    #------------------------------------------------------------------------------
    
    for(j in 1:length(stages)){ # Early marine and marine rearing
      
      # If there are timing data
      if(length(which(timing$cuid == cuid[i] & timing$stage == stages[j])) == 0){
        warning(paste0("No timing data for cuid ", cu_list$cuid[i], " ", stages[j]))
        
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
        
        # Plot flow
        if(3 == 2){
          
          
          
        }
        
        # Plot temperature
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
        
      } # end if timing data
      # print(".")
      
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