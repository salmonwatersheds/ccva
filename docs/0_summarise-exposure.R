###############################################################################
# Code to read in freshwater and marine output and write summary files for
# figures
#
# Steph Peacock
# 
###############################################################################

library(dplyr)
library(abind)

#------------------------------------------------------------------------------
# List of GCMs being considered
#------------------------------------------------------------------------------
gcms <- c("ACCESS1-0", "CanESM2", "CCSM4", "CNRM-CM5", "HadGEM2-ES", "MPI-ESM-LR")

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
cuid <- cu_list$cuid[order(cu_list$species_pooled)] # Keep same order as in output
n.CUs <- length(cuid)

# Order CU list to match incl.stages order (CK, CM, CO, PKE, PKO, SEL, SER, SH)
cu_list <- cu_list[match(cuid, cu_list$cuid), ]

#------------------------------------------------------------------------------
# Stages
#------------------------------------------------------------------------------
stages <- c("eggs_alevin", "fw_rearing", "early_marine", "marine_rearing", "adult_migration", "spawning") 

#------------------------------------------------------------------------------
# Timing
#------------------------------------------------------------------------------

# Use timing data that shortens migration window for adult lake-type sockeye to the time taken to migrate to the rearing lake, assuming that salmon holding in the lake can moderate temperatures by moving to deeper, cooler water in the lake.
timing <- read.csv("data/timing/timing-fraser.csv")
timing_days <- read.csv("data/timing/timing-fraser_days.csv")

# # Draft V1
# timing <- read.csv("data/timing/archive/timing-fraser_2024-02-19.csv")
# timing_days <- read.csv("data/timing/archive/timing-fraser_days_2024-02-19.csv")


###############################################################################
# Output 1: Combining FW and marine output on median, min, max across grid cells 
#           for each model
###############################################################################

# Change to dataframe
# CU (60), stage (6), variable (2 each stage), model (6), rcp (2), period (4) 
out.dat <- data.frame(
  species = rep(cu_list$species_abbr, each = 6*2*6*2*4),
  cu_name = rep(cu_list$cu_name_pse, each = 6*2*6*2*4),
  cuid = rep(cu_list$cuid, each = 6*2*6*2*4),
  stage = rep(rep(stages, each = 2*6*2*4), n.CUs),
  variable = rep(rep(c(1,2), each = 6*2*4), n.CUs*6),
  model = rep(rep(gcms, each = 2*4), n.CUs*6*2),
  rcp = rep(rep(c("rcp45", "rcp85"), each = 4), n.CUs*6*2*6),
  period = rep(c("hist", "early", "mid", "late"), n.CUs*6*2*6*2),
  median = NA,
  min = NA,
  max = NA
)

for(r in 1:2){
  rcp <- c("rcp45", "rcp85")[r]
  
  fw_out <- readRDS(paste0("freshwater/output/freshwater_output_fraser_", rcp, "_2024-12-19_shortSELmig.rds"))
  mar_out <- readRDS(paste0("marine/output/marine_output_fraser_", rcp, "_2024-02-06.rds"))
  
  # dim(fw_out)
  # dim(mar_out)
  # Check that order matches 
  if(sum(as.numeric(dimnames(fw_out)[[2]]) - cuid, as.numeric(dimnames(mar_out)[[2]]) - cuid) > 0){
    stop("CU orders differ!!")
  }
  
  # Paste fw and mar output along stage dimension
  out <- abind(fw_out, mar_out, along = 3)
  # dimnames(out)[[3]]
  
  # Rearrange stages in order of lifecycle
  out <- out[, , match(stages, dimnames(out)[[3]]), , ,]
  # Check: cbind(stages, dimnames(out)[[3]])
  
  # Dimensions of out
  # 1) model (6)
  # 2) CU (60)
  # 3) stage (6)
  # 4) stream temp/SST, stream flow/SSS (2)
  # 5) period (4)
  # 6) median, lower, upper through space (3)
  
  for(i in 1:n.CUs){ # For each CU
    for(j in 1:length(stages)){ # For each stage
      for(v in 1:2){ # for each climate variable (within that stage)
        for(m in 1:length(gcms)){ # For each GCM 
          ind <- which(out.dat$cuid == cuid[i] & out.dat$stage == stages[j] & out.dat$variable == v & out.dat$model == gcms[m] & out.dat$rcp == rcp)
          
          #Check
          if(length(ind) != 4){
            stop("# of rows of out.dat not equal four periods!!")
          }
          
          out.dat[ind, c("median", "min", "max")] <- out[m, i, j, v, , ]
          
          rm(ind)
        } # end m
      } # end v
    } # end j
  } #end i
} # end r

# Checks
sum(is.na(out.dat$median)) # Should be zero; no NAs
length(which(out.dat$median > out.dat$max)) # Should be zero
length(which(out.dat$median < out.dat$min)) # Should be zero

# Rename variables depending on stage
out.dat$variable[which(out.dat$stage %in% c("early_marine", "marine_rearing") & out.dat$variable == 1)] <- "SST"
out.dat$variable[which(out.dat$stage %in% c("early_marine", "marine_rearing") & out.dat$variable == 2)] <- "SSS"
out.dat$variable[which(out.dat$stage %in% c("eggs_alevin", "fw_rearing", "adult_migration", "spawning") & out.dat$variable == 1)] <- "stream_temp"
out.dat$variable[which(out.dat$stage %in% c("eggs_alevin", "fw_rearing", "adult_migration", "spawning") & out.dat$variable == 2)] <- "stream_flow"

unique(out.dat$variable)


# write.csv(out.dat, paste0("docs/ccva-output_dataframe_fraser_", Sys.Date(), "_shortSELmig.csv"), row.names = FALSE)
write.csv(out.dat, paste0("output/archive/exposure-raw_", Sys.Date(), "_shortSELmig.csv"), row.names = FALSE)
write.csv(out.dat, paste0("output/ignore/exposure-raw.csv"), row.names = FALSE)

###############################################################################
# Output 2: Exposure within and among stages
###############################################################################


# exposure_stage dims: cuid (60), stage (6), model (6), rcp (2), period (4)
exposure_stage <- data.frame(
  cuid = rep(cu_list$cuid, each = 6*6*2*4),
  stage = rep(rep(stages, each = 6*2*4), n.CUs),
  model = rep(rep(gcms, each = 2*4), n.CUs*6),
  rcp = rep(rep(c("rcp45", "rcp85"), each = 4), n.CUs*6*6),
  period = rep(c("hist", "early", "mid", "late"), n.CUs*6*6*2),
  p_j = NA,
  n_j = NA
)

# exposure dims: cuid (60), model (6), rcp (2), period (4)
exposure <- data.frame(
  cuid = rep(cu_list$cuid, each = 6*2*4),
  model = rep(rep(gcms, each = 2*4), n.CUs),
  rcp = rep(rep(c("rcp45", "rcp85"), each = 4), n.CUs*6),
  period = rep(c("hist", "early", "mid", "late"), n.CUs*6*2),
  exp_weighted = NA,
  exp_simple = NA
)

for(i in 1:n.CUs){ # for each CU
  for(m in 1:length(gcms)){ # for each model
    for(r in 1:2){ # for each emissions scenario Rcp
      for(p in 1:4){ # for each period
        
        # Average across variables within a life stage
        out.dat.ind <- which(out.dat$cuid == cuid[i] & out.dat$model == gcms[m] & out.dat$rcp == c("rcp45", "rcp85")[r] & out.dat$period == c("hist", "early", "mid", "late")[p])
        
        #Check 
        if(length(out.dat.ind) != 12){
          stop("length of out.dat.ind should equal 12 (6 stages x 2 variables")
        }
        
        # Proportion (p_j) for stage j above/below threshold averaged between two climate variables
        # Note, tapply will order output alphabetically, need to re-arrange to name
        p_j0 <- tapply(out.dat$median[out.dat.ind], out.dat$stage[out.dat.ind], mean)
        p_j <- p_j0[match(stages, names(p_j0))]
        
        # for pink and chum, set p_j to zero
        
        # Put into exposure_stage dataframe
        exposure_stage.ind <- which(exposure_stage$cuid == cuid[i] & exposure_stage$model == gcms[m] & exposure_stage$rcp == c("rcp45", "rcp85")[r] & exposure_stage$period == c("hist", "early", "mid", "late")[p])
        exposure_stage$p_j[exposure_stage.ind] <- p_j
        
        # Simple exposure is average across all p_j
        exposure.ind <- which(exposure$cuid == cuid[i] & exposure$model == gcms[m] & exposure$rcp == c("rcp45", "rcp85")[r] & exposure$period == c("hist", "early", "mid", "late")[p])
        
        exposure$exp_simple[exposure.ind] <- mean(p_j)
        
        # Weighted by time in each stage
        n_j <- timing_days$duration_total[timing_days$cuid == cuid[i]]
        names(n_j) <- timing_days$stage[timing_days$cuid == cuid[i]]
        
        exposure_stage$n_j[exposure_stage.ind] <- n_j
        
        exposure$exp_weighted[exposure.ind] <- sum(p_j * n_j)/sum(n_j)
        
        rm(p_j, p_j0, n_j, exposure.ind, exposure_stage.ind, out.dat.ind)
      } #end p
    } # end r
  } # end m
} # end i

write.csv(exposure_stage, paste0("output/archive/exposure_stages_", Sys.Date(), "_shortSELmig.csv"))
write.csv(exposure_stage, "output/ignore/exposure_stages.csv")

write.csv(exposure, paste0("output/archive/exposure_overall_", Sys.Date(), "_shortSELmig.csv"))
write.csv(exposure, "output/exposure_overall.csv")

# ###############################################################################
# # Exposure spatial
# # Summarizing across models for each grid cell
# # Combining freshwater and marine
# ###############################################################################
# 
# # Load grid cell ids for each CU and life stage
# incl.stages_fw <- readRDS("freshwater/output/PCIC_incl.rds")
# 
# # Load marine grid points
# incl.stages_em <- read.csv("marine/output/incl_NOAA_EM.csv")
# incl.stages_mr <- read.csv("marine/output/incl_NOAA_marRear.csv")
# 
# # Create list to store full spatial output for plotting maps
# # Later will calculate median across all GCMs
# exposure_spat <- list(); length(exposure_spat) <- n.CUs*length(stages)*2*4
# dim(fw_spat) <- c(n.CUs, length(stages))
# dimnames(fw_spat) <- list(cuid, stages)
# 
# for(i in 1:n.CUs){
#   for(j in 1:length(stages)){
#     fw_spat[[i, j]] <- array(data = NA,
#                              dim = c(n.models, 2, 4, length(incl.stages[[i,j]])),
#                              dimnames = list(
#                                gcms$modelName,
#                                c("temp", "flow"),
#                                c("hist", "early", "mid", "late"),
#                                incl.stages[[i,j]])
#     )
#   }
# }
# 
# fw_spat <- readRDS(paste0("freshwater/output/freshwater_spat_fraser_rcp", rcp, "_2024-06-26.rds"))
# 
# # write.csv(out.spat, paste0("output/exposure-spatial_fraser_", Sys.Date(), "_shortSELmig.csv"), row.names = FALSE)
# 
# 
# ###############################################################################
