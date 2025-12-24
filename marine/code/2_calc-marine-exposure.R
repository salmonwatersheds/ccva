###############################################################################
# Code to summarize gridded climate change projections for sea surface temp.
# and salinity from CMIP5 GCMs and calculate projected changes in exposure to
# temperature and salinity outside of the optimal range for each species 
# (for temperature) and historical norms in each grid cell (for salinity).
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

# Create pooled species field (for matching to temp max)
cu_list$species_pooled <- cu_list$species_name
cu_list$species_pooled[cu_list$species_name %in% c("Lake sockeye", "River sockeye")] <- "Sockeye"
cu_list$species_pooled[cu_list$species_name %in% c("Pink (odd)", "Pink (even)")] <- "Pink"

# Create vector of cuid
cuid <- cu_list$cuid[order(cu_list$species_pooled)] # Keep same order as in freshwater
n.CUs <- length(cuid)

# Order CU list to match incl.stages order (CK, CM, CO, PKE, PKO, SEL, SER, SH)
cu_list <- cu_list[match(cuid, cu_list$cuid), ]


#------------------------------------------------------------------------------
# Timing
#------------------------------------------------------------------------------

timing <- read.csv("data/timing/timing-fraser.csv") %>%
  filter(stage %in% c("early_marine", "marine_rearing"))

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
    "CAN-ESM2",
    "CCSM4",
    "CNRM-CM5",
    "HADGEM2-ES",
    "MPI-ESM-LR")
)

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
    model = gcms$modelName[m], # Which GCM?
    rcp = rcp, # Which emissions scenario? Options 45 or 85
    spat_id = spat_overall 
  )
  
  SSS <- loadNOAA(
    variable = "SSS", # Which variable to load? One of "PP", "pH", "SSS", "SST"
    model = gcms$modelName[m], # Which GCM? 
    rcp = rcp, # Which emissions scenario? Options 45 or 85
    spat_id = spat_overall 
  )
  
  
  # Check dimnames: this is how output will be organized
  grid.ref <- as.numeric(dimnames(SST)[[1]])
  
  #------------------------------------------------------------------------------
  # Loop through each CU
  #------------------------------------------------------------------------------
  
  for(i in 1:n.CUs){
    
    #------------------------------------------------------------------------------
    # Loop through each life stage
    #------------------------------------------------------------------------------
    
    for(j in 1:length(stages)){ # Early marine and marine rearing
      
        # Extract timing and covert to day-of-year
        timing.ij <- timing[which(timing$cuid == cuid[i] & timing$stage == stages[j])[1], c("start", "duration")]
        
        #-----------------------------------------------------------------------------
        # Identify grid cells for given life stage
        #-----------------------------------------------------------------------------
        
        # See freshwater-grid.R for code that identifies which grid cells should be used
        # for each life stage.
        
        incl <- list(spat_EM, spat_MR)[[j]]
        
        # # Are there missing data?
        # nas <- which(is.na(SST[match(incl, grid.ref), ]), arr.ind = TRUE)
        # length((unique(nas[, 1]))) # Slightly more NAs for SSS than SST when m = 1
        # 
        # Dimensions
        n.grid <- length(incl)
        
        #-----------------------------------------------------------------------------
        # Subset stream temp relevant to CU and stage
        #-----------------------------------------------------------------------------
        T.ij <- S.ij <- array(NA,
                                dim = c(n.grid, 4, timing.ij$duration, 30),
                                dimnames = list(incl, c("hist", "early", "mid", "late"), NULL, NULL))
        # Note: for stages that span >365 days, some temps may be counted twice, because essentially two cohorts occupy that space at the same time.
        
        for(p in 1:4){ # for each period
          for(y in 1:30){ # for each year in the period
            
            start.ind <- which(years == period.years[y, p] & months == timing.ij$start)
            
            # For the last year in late century, can't extend if start + duration > 365
            if((start.ind + timing.ij$duration - 1) > dim(SST)[2]){
              
              # n.available <- length(start.ind:(dim(SST)[2]))
              # T.ij[, p, 1:n.available, y] <- SST[match(incl, grid.ref), start.ind:(dim(SST)[2])]
              # DO NOTHING; leave all as NA otherwise biases period averages?
              
            } else {
              
              T.ij[, p, , y] <- SST[match(incl, grid.ref), start.ind:(start.ind + timing.ij$duration - 1)]
              
            }
            
            # For the last year in late century, can't extend if start + duration > 365
            # Separate check for SSS because dims aren't always consistent for HADGEM2-ES
            if((start.ind + timing.ij$duration - 1) > dim(SSS)[2]){
              
              # n.available <- length(start.ind:(dim(SSS)[2]))
              # S.ij[, p, 1:n.available, y] <- SSS[match(incl, grid.ref), start.ind:(dim(SSS)[2])]
              # DO NOTHING; leave all as NA otherwise biases period averages?
              
              
            } else {
              
             S.ij[, p, , y] <- SSS[match(incl, grid.ref), start.ind:(start.ind + timing.ij$duration - 1)]
            }
            
          }}
        
  
        # Plot sST
        if(3 == 2){
          xx <- 1 # which grid cell
          
          plot(st_geometry(st_buffer(em_zone, dist = 100*1000)), col = NA, border = NA)
          plot(st_geometry(em_zone), col = cols[5], border = NA, add = TRUE)
          plot(st_geometry(grid_points[match(incl[xx], grid_points$id), ]), add = TRUE, 1, pch = 19)
          
          start.ind <- which(years == period.years[y, p] & months == timing.ij$start)
          
          # Plot 5 years
          par(mar = c(4,4,4,1), bg = "white")
          plot(date[years %in% period.years[1:5, 4]], SST[match(incl[xx], grid.ref), which(years %in% period.years[1:5, 4])], "l", xlab = "Year", ylab = "SST (˚C)", las = 1, bty = "l", xaxs ="i", col = cols[1], lwd = 1.5)
          for(p in 1:3){
            lines(date[years %in% period.years[1:5, 4]], SST[match(incl[xx], grid.ref), which(years %in% period.years[1:5, p])], lwd = 1.5, col = cols[c(5,3,4)[p]])
          }
          abline(h = tmax[i], lty = 2)
          
          for(p in 1:4){
            for(y in 1:5){
            start.ind <- which(years %in% period.years[y, 4] & months == timing.ij$start)
            points(date[start.ind:(start.ind + timing.ij$duration - 1)], T.ij[xx, p, , y], pch = 19, cex =0.8, col = cols[c(5,3,4,1)[p]])
          }}
          
      # Plot salinity
          xxx <- 50578 #middle of Pacific
          # xxx <- incl[xx]
          plot(date[years %in% period.years[10:14, 3]], SSS[match(xxx, grid.ref), which(years %in% period.years[10:14, 3])], "o", xlab = "Year", ylab = "Salinity", las = 1, bty = "l", xaxs ="i", col = cols[3], lwd = 1.5, pch = 19)
          lines(date[years %in% period.years[10:14, 3]], SSS[match(xxx, grid.ref), which(years %in% period.years[10:14, 1])], "o", col = cols[5], pch = 19, xpd = NA)
          abline(h = quantile(SSS[match(xxx, grid.ref), which(years %in% period.years[, 1])], 0.025), lty = 2)
          
          hist(SSS[match(xxx, grid.ref), which(years %in% period.years[, 1])], col = cols[5], main = "Historical", xlab = "Salinity (ppt)", yaxs = "i", border = "white", breaks = seq(31.5, 33, 0.1))
          hist(SSS[match(xxx, grid.ref), which(years %in% period.years[, 3])], col = cols[3], add = TRUE, breaks = seq(31.5, 33, 0.1), xpd = NA)
          abline(v = quantile(SSS[match(xxx, grid.ref), which(years %in% period.years[, 1])], 0.025), lty = 2, lwd = 2)
          
        } # end if 3 == 2 plot
        
        # Collapse Ts.ij years and days for easy
        dim(T.ij) <- c(n.grid, 4, 30*timing.ij$duration) 
        dim(S.ij) <- c(n.grid, 4, 30*timing.ij$duration) 
        
        #-----------------------------------------------------------------------------
        # Sea surface temperature 
        #-----------------------------------------------------------------------------
        
        #----------------------------------------------------------------------
        # Calculate months above threshold and % DAT
        #----------------------------------------------------------------------
        T.mat <- T.perc <- array(NA, 
                                   dim = c(4, n.grid),
                                   dimnames = list(c("hist", "early", "mid", 'late'), NULL)
        )
        
        for(p in 1:4){ # For each period
          
          
          T.mat[p, ] <- apply(T.ij[, p, ] > matrix(rep(tmax[i], 30*timing.ij$duration), n.grid, 30*timing.ij$duration), 1, sum, na.rm = TRUE)
          
          
          # Percentage of days;
          # Need to account for late-century period not being able to extend into next year for 2099 if timing.ij$start + timing.ij$n.days > 365
          # Count number of months that are not NA
          n_months_nna <-  apply(T.ij[, p, ], 1, function(x){sum(!is.na(x))})
          T.perc[p, ] <- T.mat[p, ] / n_months_nna
          T.perc[p, which(n_months_nna == 0)] <- NA # Don't keep -Inf
          rm(n_months_nna)
        } # end p
        
        # Store full spatial output
        mar_spat[[i,j]][m, "SST", , ] <- T.perc
        
        # store summary output
        for(p in 1:4){ # for each period
          # Store output in end array
          mar_output[m, i, j, "SST", p, ] <- c(
            median(T.perc[p, ], na.rm = TRUE),
            range(T.perc[p, ], na.rm = TRUE))
          
        } # end p
        
        #-----------------------------------------------------------------------------
        # Salinity
        #-----------------------------------------------------------------------------
        
        # Calculate 2.5th percentile over historical period for each included grid cell
        # (Use all DOY, not Qs.ij for life stage)
        hist_SSS <- apply(SSS[match(incl, grid.ref), which(period == "hist")], 1, quantile, 0.025, na.rm = TRUE)
        
        # # What does salinity look like, historically
        # # Early marine xx <- 1182; marine rearing xx <- 2019
        # plot(date[which(period == "hist")], SSS[xx, which(period == "hist")], "l", ylim= c(30.5, 33.5))
        # lines(date[which(period == "hist")], SSS[xx, which(period == "early")], col = cols[5]) 
        # lines(date[which(period == "hist")], SSS[xx, which(period == "mid")], col = cols[3]) 
        # lines(date[which(period == "hist")], SSS[xx, which(period == "late")], col = cols[1]) 
        # abline(h = quantile(SSS[xx, which(period == "hist")], 0.025, na.rm = TRUE))
        
        # Calculate days below 2.5th percentile of historical 
        S.dbt <- S.perc <- array(NA, 
                                   dim = c(4, n.grid),
                                   dimnames = list(c("hist", "early", "mid", 'late'), 
                                                   NULL)
        )
        
        for(p in 1:4){
          S.dbt[p, ] <- apply(S.ij[, p, ] < matrix(rep(hist_SSS, 30*timing.ij$duration), n.grid, 30*timing.ij$duration), 1, sum, na.rm = TRUE)
          
          # Percentage of days;
          # Need to account for late-century period not being able to extend into next year for 2099 if timing.ij$start + timing.ij$n.days > 365
          n_months_nna <-  apply(S.ij[, p, ], 1, function(x){sum(!is.na(x))})
          
          S.perc[p, ] <- S.dbt[p, ] / n_months_nna
          S.perc[p, which(n_months_nna == 0)] <- NA # Don't keep -Inf
          rm(n_months_nna)
          
        } # end p
        
        # # Check
        # plot(c(1:4), S.perc[, 14])
        # 
        
        # Store full spatial output
        mar_spat[[i,j]][m, "SSS", , ] <- S.perc
        
        for(p in 1:4){ # for each period
          # Store output in end array
          mar_output[m, i, j, "SSS", p, ] <- c(
            median(S.perc[p, ], na.rm = TRUE),
            range(S.perc[p, ], na.rm = TRUE))
        }
        
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

saveRDS(mar_spat, file = paste0("marine/output/marine_spat_fraser_rcp", rcp, "_", Sys.Date(), ".rds"))
saveRDS(mar_output, file = paste0("marine/output/marine_output_fraser_rcp", rcp, "_",  Sys.Date(), ".rds"))
