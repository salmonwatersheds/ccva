###############################################################################
# Code to load PCIC model output, trim to spatial and temporal domain of interest,
# calculate weekly rolling max or average, and save output RDS file
# for use in calc-freshwater-exposure.R
# 
# Contact: Stephanie Peacock <speacock@psf.ca>
# Date: Aug 31, 2023
###############################################################################

source('code/load-PCIC-Fraser.R')

###############################################################################
# Possible Global Climate Models (GCMs)
###############################################################################

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
# Loop through each GCM and save output
###############################################################################

for(m in 2:n.models){ # for each model
  
  #----------------------------------------------------------------------------
  # Stream temperature
  #----------------------------------------------------------------------------
  # Load model output; this takes a minute
  Ts <- loadPCIC(
    variable = "waterTemp", # Which variable to load? One of waterTemp or discharge
    model = gcms$modelName[m] # Which GCM?
  )
  
  # Change Ts from degrees Kelvin to Celsius
  if(min(Ts[, 1000], na.rm = TRUE) > 200){
    Ts <- Ts - 273.15
  }
  
  # Calculate 7-day rolling mean or max (takes 2.3 mins)
  # start.time <- Sys.time()
  Ts.weeklyMax <- t(apply(Ts, 1, weeklyStat, stat = "max"))
  # print(Sys.time() - start.time)
  
  
  #----------------------------------------------------------------------------
  # Stream flow
  #----------------------------------------------------------------------------
  
  Qs <- loadPCIC(
    variable = "discharge",
    model = gcms$modelName[m]
  )
  
  Qs.weeklyMean <- t(apply(Qs, 1, weeklyStat, stat = "mean"))
  
  #----------------------------------------------------------------------------
  # Combine into single object and trim extra week used for rolling mean calc
  #----------------------------------------------------------------------------
  
  dates <- dimnames(Ts.weeklyMax)[[2]]
  var <- list(
    Ts = Ts.weeklyMax[, which(dates == "1970-01-01"):which(dates == "2099-12-31")],
    Qs = Qs.weeklyMean[, which(dates == "1970-01-01"):which(dates == "2099-12-31")]
  )
  
  #----------------------------------------------------------------------------
  # Save output
  #----------------------------------------------------------------------------
  
  saveRDS(var, paste0("data/processed-data/PCIC_", gcms$modelName[m], "_processed.rds"))
  
  print(paste0("End ", gcms$modelName[m]))
} #end m
  