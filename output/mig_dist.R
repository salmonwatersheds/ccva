# Are Chinook exposure ranked interior to exterior?

# median overall exposure
me <- exp_summary0[match(cu_order, as.numeric(names(me)))[1:19],1]

cu_center <- st_centroid(cu_boundary_ck)

cu_boundary_ck <- cu_boundary %>% filter(species == "Chinook")
cu_boundary_ck <- cu_boundary_ck[match(cu_order, cu_boundary_ck$CUID)[1:19],]

plot(st_geometry(cu_boundary_ck), col = col_palette_fun(me))

# Migration distance across all CUs
dat_root <- "freshwater/data/spatial/"
streams <- st_read(paste0(dat_root, "bcfishpass.gdb"), layer = "streams") %>%
  st_transform(crs = 4269)

# Create lookup for spawning and rearing fields in the geodatabase
spp_lookup <- data.frame(
  species_pooled = sort(unique(cu_list$species_pooled)),
  streams_code = c("ch", "cm", "co", "pk", "sk", "st")
)


ck_cuid <- cu_order[1:19]
mig_dist <- rep(NA, 60)
ind <- 0
for(s in 5:6){
  
  # Load migration paths for that species
  ind.mig <- st_read(paste0(dat_root, "bcfishpass.gdb"), layer = paste0("cu_migrationpaths_", spp_lookup$streams_code[s]))
  
  # Extract cuids for the selected species
  cuid.s <- cu_list$pooledcuid[which(cu_list$species_pooled == spp_lookup$species_pooled[s])]
  
  for(i in 1:length(cuid.s)){
    
    if(cuid.s[i] %in% unique(ind.mig$cuid)){
      
      streams.i <- streams[streams$linear_feature_id %in% ind.mig$linear_feature_id[ind.mig$cuid == cuid.s[i]], ]
      plot(st_geometry(streams.i))
      plot(st_geometry(cu_boundary[cu_boundary$CUID == cuid.s[i], ]), add = TRUE, col = NA, border = 2)
      
      
    mig_dist[which(cu_list$cuid == cuid.s[i])] <- c( %>% select(length_metre) %>% summarise(sum(length_metre)*10^-3) %>% as.numeric())[1]
    }
    
  } #end i
}

plot(mig_dist, exp_summary0[, 1], pch = as.numeric(as.factor(cu_list$species_pooled)), col = species_cols[cu_list$species_pooled], cex = 1.5, lwd = 2)

fit <- lm(exp_summary0[, 1] ~ cu_list$species_pooled + mig_dist)
summary(fit)


# Chinook exposure summer vs fall
ck_order <- data.frame(
  cuid = c(
  c(cu_list$cuid[which(cu_list$species_abbr == "CK" & grepl("Spring", cu_list$cu_name_pse))],
cu_list$cuid[which(cu_list$species_abbr == "CK" & grepl("Summer", cu_list$cu_name_pse))],
cu_list$cuid[which(cu_list$species_abbr == "CK" & grepl("Fall", cu_list$cu_name_pse))]))
)

ck_order$mig_timing <- c(
  rep("Spring", length(which(cu_list$species_abbr == "CK" & grepl("Spring", cu_list$cu_name_pse)))),
  rep("Summer", length(which(cu_list$species_abbr == "CK" & grepl("Spring", cu_list$cu_name_pse)))),
  rep("Fall", length(which(cu_list$species_abbr == "CK" & grepl("Spring", cu_list$cu_name_pse))))
)

ck_order$rear_type <- NA
ck_order$rear_type[which((grepl("(4-2)", cu_list$cu_name_pse[match(ck_order$cuid, cu_list$cuid)]) | grepl("(5-2)", cu_list$cu_name_pse[match(ck_order$cuid, cu_list$cuid)])))] <- "Stream"
ck_order$rear_type[grep("(4-1)", cu_list$cu_name_pse[match(ck_order$cuid, cu_list$cuid)])] <- "Ocean"


ck_order$cu_name <- cu_list$cu_name_pse[match(ck_order$cuid, cu_list$cuid)]
ck_order$overall_exposure <- exp_summary0[match(ck_order$cuid, cu_list$cuid), 1]

summary(lm(overall_exposure ~ mig_timing * rear_type, data = ck_order))

z <- stagevar_exposure %>% filter(species == "CK", stage %in% c("adult_migration", "spawning"), variable == "stream_temp")

z$mig_timing <- ck_order$mig_timing[match(z$cuid, ck_order$cuid)]
z$rear_type <- ck_order$rear_type[match(z$cuid, ck_order$cuid)]

summary(lm(exp_median ~ mig_timing*stage + rear_type* stage, data = z))

par(mfrow = c(1,2))
plot(factor(z$mig_timing[z$stage == "adult_migration"], levels = c("Spring", "Summer", "Fall")), z$exp_median[z$stage == "adult_migration"], pch = 2, col = 4, main = "Adult migration", ylab = "Exposure", xlab = "Migration timing")

plot(factor(z$mig_timing[z$stage == "spawning"], levels = c("Spring", "Summer", "Fall")), z$exp_median[z$stage == "spawning"], pch = 2, col = 3, main = "Spawning", ylab = "Exposure", xlab = "Migration timing")
