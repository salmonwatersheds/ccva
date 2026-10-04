# Prep the spatial data for showing in the app

spat_fw <- readRDS("freshwater/output/freshwater_spat_fraser_rcp45_2024-12-19_shortSELmig.rds")

# dim 60 x 4
# CU by life stage
# Wihtin each element, dim 6, 2, 4, 9
# (6) GCM
# (2) variables
# (4) periods
# grid cells (variable)

spat_mar <- readRDS("marine/output/marine_spat_fraser_rcp45_2024-02-06.rds")
# 60 x 2


# We want an output for each CU and life stage and variable, mid-century, averaged across models.
# List 60 x 6
# 
stages <- c("eggs_alevin", "fw_rearing", "early_marine", "marine_rearing", "adult_migration", 'spawning') 
spat <- list(); length(spat) <- 60*6; dim(spat) <- c(60, 6)
dimnames(spat) <- list(
  dimnames(spat_mar)[[1]],
  stages
)

for(i in 1:60){
  for(j in 1:6){
    if(stages[j] %in% dimnames(spat_fw)[[2]]){
      spat.ij <- spat_fw[[i, match(stages[j], dimnames(spat_fw)[[2]])]]
    } else {
      spat.ij <- spat_mar[[i, match(stages[j], dimnames(spat_mar)[[2]])]]
    }
    
    spat[[i,j]] <- array(NA, dim = c(2, dim(spat.ij)[4]))
    for(v in 1:2){ # for each variable
      # Take median across models
      if(dim(spat.ij)[4] == 1) {
        spat[[i,j]][v, ] <- median(spat.ij[ , v, 3, ])
      } else {
        spat[[i,j]][v, ] <- apply(spat.ij[ , v, 3, ], 2, median, na.rm = TRUE)
      }
    } # end v
  } # end j
} # end i

# Check
# spat[[47, 5]][1,]
saveRDS(spat, "shiny/data/spat_combined.rds")


#------------------------------------------------------------------------------
# Overall exposure
#------------------------------------------------------------------------------
exposure <- read.csv("output/exposure_overall.csv") # Summarised across stages

ind <- which(exposure$rcp == 'rcp45' & exposure$period == "mid")
var <- which(names(exposure) == c("exp_simple"))
exp_summary0 <- tapply(exposure[ind, var], exposure$cuid[ind], median)


# re-order by cu_list cuids
dum.cuid <- as.numeric(names(exp_summary0))
exp_summary <- exp_summary0[match(cu_list$cuid, dum.cuid)]

write.csv(exp_summary, "shiny/data/overall_exposure_rcp45_mid.csv", row.names = FALSE)

#------------------------------------------------------------------------------
# Stage exposure
#------------------------------------------------------------------------------

exposure_stage <- read.csv("output/ignore/exposure_stages.csv") %>%
  filter(rcp == "rcp45", period == "mid") %>%
  mutate(cu_stage = paste(cuid, stage)) %>%
  group_by(cu_stage) %>%
  summarise(cuid = unique(cuid),
            stage = unique(stage),
            p_j_median = median(p_j)) %>%
  select(cuid, stage, p_j_median)

# Change pink and chum to NA

write.csv(exposure_stage, "shiny/data/stage_exposure_rcp45_mid.csv", row.names = FALSE)

#------------------------------------------------------------------------------
# Stage and variable exposure
#------------------------------------------------------------------------------
exposure <- read.csv("output/archive/exposure-raw_2024-12-20_shortSELmig.csv") %>%
  filter(rcp == "rcp45", period == "mid") %>%
  mutate(cu_stage_var = paste(cuid, stage, variable)) %>%
  group_by(cu_stage_var) %>%
  summarise(
    species = unique(species),
    cuid = unique(cuid),
    cu_name = unique(cu_name),
    stage = unique(stage),
    variable = unique(variable),
    exp_median = median(median),
    exp_lowestGCM = min(median),
    exp_highestGCM = max(median)) %>%
  select(species, cuid, cu_name, stage, variable, exp_median, exp_lowestGCM, exp_highestGCM)

write.csv(exposure, "output/stagevariable_exposure_rcp45_mid.csv", row.names = FALSE)
