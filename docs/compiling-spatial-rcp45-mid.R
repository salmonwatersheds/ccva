###############################################################################
# Compile spatial output for RCP 4.5 and mid-century
# across marine and freshwater for mapping in Shiny app
# Feb 19, 2024
###############################################################################

mar_spat <- readRDS("marine/output/marine_spat_fraser_rcp45_2024-02-06.rds")
fw_spat <- readRDS("freshwater/output/freshwater_spat_fraser_rcp45_2024-02-12.rds")

# dimensions:
#   n.CUs,
#  stages

dimnames(fw_spat[[1,1]])
# n.models = 6
# temp, flow variables = 2
# hist/early/mid/late = 4
# grid id = xxx (different for each CU)

# organize into combined list with median across all models
stages <- c("eggs_alevin", "fw_rearing", "early_marine", "marine_rearing", "adult_migration", "spawning")

all_spat <- list(); length(all_spat) <- n.CUs*length(stages)
dim(all_spat) <- c(n.CUs, length(stages))
dimnames(all_spat) <- list(dimnames(fw_spat)[[1]], stages)

for(i in 1:n.CUs){
  for(j in 1:length(stages)){
    
    if(j %in% c(1,2,5,6)){ # Freshwater
      
      dum <- fw_spat[[i, match(stages[j], dimnames(fw_spat)[[2]])]][, , 3, ]
      if(length(dim(dum)) == 2){ # Iff just one grid cell
        all_spat[[i, j]] <- apply(dum, 2, median, na.rm = TRUE)
      } else {
        all_spat[[i, j]] <- array(NA, dim = c(dim(dum)[2:3]), dimnames = list(dimnames(dum)[[2]], dimnames(dum)[[3]]))
        for(k in 1:2){
          all_spat[[i, j]][k, ] <- apply(dum[, k, ], 2, median, na.rm = TRUE)
        }
      }
      rm(dum)
      
    } else { # Marine
      
      dum <- mar_spat[[i, match(stages[j], dimnames(mar_spat)[[2]])]][, , 3, ]
      all_spat[[i, j]] <- array(NA, dim = c(dim(dum)[2:3]), dimnames = list(dimnames(dum)[[2]], dimnames(dum)[[3]]))
      for(k in 1:2){
        all_spat[[i, j]][k, ] <- apply(dum[, k, ], 2, median, na.rm = TRUE)
      }
      rm(dum)
    }
  }}

saveRDS(all_spat, file = "docs/data/all_spat.RDS")
