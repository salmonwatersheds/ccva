###############################################################################
#
# 1a_spatial-distribution.R
#
# This code reads in a number of spatial datasets related to salmon distribution
# in the **Fraser region** and determines, for each Conservation Unit (CU),  
# which PCIC grid cells should be used in assessments of climate change exposure  
# for freshwater life stages.
#
# Contact: Steph Peacock (speacock@psf.ca)
# Date: Jan 3, 2024
###############################################################################

library(dplyr)
library(sf)
library(wesanderson)

# Colour palette for plotting
cols <- wes_palette("Darjeeling1")

# Set root for spatial datasets
dat_root <- "freshwater/data/spatial/"

# Set root for X Drive (user dependent; assuming ccva repo is in X Drive/1_PROJECTS)
XDrive_root <- paste(strsplit(getwd(), "/")[[1]][1:6], collapse = "/")

###############################################################################
# Read in relevant spatial datasets
###############################################################################

#------------------------------------------------------------------------------
# Up-to-date CU list (taken from database)
#------------------------------------------------------------------------------

cu_list <- read.csv("data/conservationunits_decoder.csv") %>%
	subset(region == "Fraser") 

# Create variable for pooled species
cu_list$species_pooled <- cu_list$species_name
cu_list$species_pooled[cu_list$species_name %in% c("Pink (odd)", "Pink (even)")] <- "Pink"
cu_list$species_pooled[cu_list$species_name %in% c("Lake sockeye", "River sockeye")] <- "Sockeye"

# Remove CUs that are considered Extinct
cu_list <- cu_list[which(cu_list$cu_type != "Extinct" | is.na(cu_list$cu_type)), ]

#------------------------------------------------------------------------------
# bcfishpass: spawning and rearing stream distributions for each species
# Nov 23, 2023: NOTE: these distributions are *usable* habitat and do not account
# for habitats that salmon may transit through.
#
# Updated to use version of database from https://www.hillcrestgeo.ca/outgoing/forPSF/
# Dated 2023-Dec-08 09:11 
#------------------------------------------------------------------------------

streams <- st_read(paste0(dat_root, "bcfishpass.gdb"), layer = "streams") %>%
	st_transform(crs = 4269)

# Create lookup for spawning and rearing fields in the geodatabase
spp_lookup <- data.frame(
	species_pooled = sort(unique(cu_list$species_pooled)),
	streams_code = c("ch", "cm", "co", "pk", "sk", "st")
)

# If encounter the error: 
# Error in wk_handle.wk_wkb(wkb, s2_geography_writer(oriented = oriented,  : 
# Loop 0 is not valid: Edge 2607 has duplicate vertex with edge 2625
# Then run:
sf_use_s2(FALSE)

# Load migration paths
# List layers in gdb
# library(rgdal)
# ogrListLayers(paste0(dat_root, "bcfishpass.gdb"))
# Load when looping through species

#------------------------------------------------------------------------------
# Conservation Unit boundaries for Fraser CUs (all species)
#------------------------------------------------------------------------------

cu_boundary <- st_read(paste0(XDrive_root, "/5_DATA/CUs_Master/GDB/PSF_CUs_Master.gdb")) %>%
  subset(regionname == "Fraser") %>%
  st_transform(crs = 4269)

# Are all CUs in cu_list in cu_boundary?
cu_list$cuid %in% cu_boundary$CUID # Yes

# #------------------------------------------------------------------------------
# # Migration paths from spawning habitat to ocean entry point
# # NOTE: these were based on NuSEDS spawner survey locations across all species
# # This analysis will need to be updated so that the bcfishpass rearing and spawning 
# # stream reaches are all connected by a mig path to ocean entry point
# # Missing completely for steelhead, so have to construct for that species
# # June 24, 2024 - this is deprecated now that we have better mig paths from Simon Norris
# #------------------------------------------------------------------------------
# mig_paths <- st_read(dsn = paste0(dat_root, "fw-migration/spawn_timing_migration_paths_wCU_fraser.shp")) %>% 
# 	st_transform(crs = 4269)
# 
# cu_list$cuid %in% mig_paths$cuid # Doesn't include steelhead
# 
# mig_paths_start <- data.frame(
#   csv_id = mig_paths$csv_id,
#   lat = mig_paths$lat_snap,
#   lon = mig_paths$lon_snap) %>%
#   st_as_sf(coords = c("lon", "lat"), crs = 4269)

#------------------------------------------------------------------------------
# PCIC grid points
#------------------------------------------------------------------------------

# Read in PCIC grid
grid_points0 <- read.csv("freshwater/data/processed-data/PCIC-grid-points_fraser.csv") 

# Create grid polys
n <- length(grid_points0$id)
d <- 1/16

grid_polys <- st_as_sf(data.frame(
	id = rep(grid_points0$id, each = 5),
	rep(c("SW0", "NW", "NE", "SE", "SW1"), n),
	lon = c(rep(grid_points0$lon, each = 5) + rep(c(- d/2, -d/2, d/2,  d/2, -d/2), n)),
	lat = c(rep(grid_points0$lat, each = 5) + rep(c(- d/2,  d/2, d/2, -d/2, -d/2), n))
), coords = c("lon", "lat"), crs = 4269) %>% 
	group_by(id) %>%
	summarise(geometry = st_combine(geometry)) %>%
	st_cast("POLYGON") 

# saveRDS(grid_polys, file = "docs/data/grid_polys_fw.rds")

# Convert grid points to spatial object 
grid_points <- st_as_sf(grid_points0, coords = c("lon", "lat"), crs = 4269)

#------------------------------------------------------------------------------
# Load shoreline data (for mapping only)
#------------------------------------------------------------------------------

shoreline <- st_read(paste0(dat_root, "layers/GSHHS_i_L1.shp")) %>%
	st_transform(crs = 4269)

###############################################################################
# Loop through CUs and life stages
###############################################################################

#-----------------------------------------------------------------------------
# Define four freshwater stages
#-----------------------------------------------------------------------------
stages <- c("adult_migration", "spawning", "eggs_alevin", "fw_rearing")

#-----------------------------------------------------------------------------
# Define Fraser CU cuid, ordered CK, CM, CO, PKO, SEL, SER, SH
#-----------------------------------------------------------------------------

cuid <- cu_list$pooledcuid[order(cu_list$species_abbr)]

#-----------------------------------------------------------------------------
# Loop through to extract grid cell ids
#-----------------------------------------------------------------------------

PCIC_incl <- vector(mode = 'list', 
									 length = length(cuid)*length(stages))
dim(PCIC_incl) <- c(length(cuid), length(stages))
dimnames(PCIC_incl) <- list(cuid, stages)


I <- 0 # Overall cu # counter

# For each species
for(s in 1:6){
  
  # Load migration paths for that species
  ind.mig <- st_read(paste0(dat_root, "bcfishpass.gdb"), layer = paste0("cu_migrationpaths_", spp_lookup$streams_code[s]))
	
	#-----
  # Subset spawning and rearing streams for that species (this is a bit slow)
	
  # SPAWNING
  ind.spawn <- streams %>% 
		as.data.frame %>% 
		select(paste0("model_spawning_", spp_lookup$streams_code[s]))
	
	spawning.s <- streams[which(ind.spawn == 1), ] # %>%  
	# filter(is.na(barriers_pscis_dnstr)) %>% # No PSCIS barrier downstream; new addition June 26, 2024
	#   filter(is.na(barriers_dams_dnstr)) # No dam downstream; new addition June 26, 2024
	# 
	# if(s < 6){ # If not steelhead, remove upstream of gradient or other natural barrier to salmon
	#   spawning.s <- spawning.s %>%
	#     filter(barriers_ch_cm_co_pk_sk_dnstr == "") # No gradient or other natural barrier downstream
	# } else {
	#   spawning.s <- spawning.s %>%
	#     filter(barriers_st_dnstr == "") # No gradient or other natural barrier downstream
	# }
	

	# REARING
	if(spp_lookup$species_pooled[s] %in% c("Chum", "Pink") == FALSE) { # No rearing distributions for pink and chum
	  ind.rear <- streams %>% 
	    as.data.frame %>% 
	    select(paste0("model_rearing_", spp_lookup$streams_code[s]))
	  
	  rearing.s <- streams[which(ind.rear == 1), ] #%>%  
	  #   filter(is.na(barriers_pscis_dnstr)) %>% # No PSCIS barrier downstream; new addition June 26, 2024
	  #   filter(is.na(barriers_dams_dnstr)) # No dam downstream; new addition June 26, 2024
	  # 
	  # if(s < 6){ # If not steelhead, remove upstream of gradient or other natural barrier to salmon
	  #   rearing.s <- rearing.s %>%
	  #     filter(barriers_ch_cm_co_pk_sk_dnstr == "") # No gradient or other natural barrier downstream
	  # } else {
	  #   rearing.s <- rearing.s %>%
	  #     filter(barriers_st_dnstr == "") # No gradient or other natural barrier downstream
	  # }
	}
	
	# # Plot all CUs for that species
	# plot(st_geometry(spawning.s), col = cols[2], axes = TRUE, ylim = c(50, 53), xlim = c(-124, -118), lwd = 2)
	# plot(st_geometry(which(cu_boundary$species == "Chinook", ]), col = "#00000030", border = "#00000060", add = TRUE)
	# plot(st_geometry(shoreline), add = TRUE)
	# plot(st_geometry(streams[which(ind.spawn == 1), ] ), col = paste0(cols[2], 30), lwd = 2, add = TRUE)
	# # plot(st_geometry(streams[which(ind.rear == 1), ] ), col = paste0(cols[1], 30), lwd = 1, add = TRUE)
	# # plot(st_geometry(rearing.s), col = cols[1], lwd = 1, add = TRUE)
	# plot(st_geometry(streams %>% filter(linear_feature_id %in% ind.mig$linear_feature_id) ), col = cols[3], lwd = 0.5, add = TRUE)
	
	# Extract cuids for the selected species
	cuid.s <- cu_list$pooledcuid[which(cu_list$species_pooled == spp_lookup$species_pooled[s])]
	# cu_list$cu_name_pse[which(cu_list$species_pooled == spp_lookup$species_pooled[s])]
	
	# For each CU of that species
	for(i in 1:length(cuid.s)){
		# Subset CU boundary
		cu_boundary.i <- cu_boundary[cu_boundary$CUID == cuid.s[i], ]
		
		# MIGRATION
		if(cu_list$cu_name_pse[which(cu_list$species_pooled == spp_lookup$species_pooled[s])[i]] == "Harrison-Upstream Migrating-Late"){
		# If upstream migrating, use mig path for downstream migrating
		  mig_paths.i <- streams[streams$linear_feature_id %in% ind.mig$linear_feature_id[ind.mig$cuid == 713], ]
		  } else {
		  mig_paths.i <- streams[streams$linear_feature_id %in% ind.mig$linear_feature_id[ind.mig$cuid == cuid.s[i]], ]
		  
		}
		
		intrscts_mig <- st_intersects(mig_paths.i, grid_polys, sparse = FALSE)
		incl_mig <- which(apply(intrscts_mig, 2, sum) > 0)
		
		# SPAWNING: Find grid_id that intersect freshwater *spawning*
		dum.spawn <- st_intersects(spawning.s, cu_boundary.i, sparse = FALSE)
		spawning.i <- spawning.s[which(dum.spawn == TRUE),]
		
		intrscts_spawn <- st_intersects(spawning.i, grid_polys, sparse = FALSE)
		incl_spawn <- which(apply(intrscts_spawn, 2, sum) > 0)
		
		# REARING: Find grid_id that intersect freshwater *rearing* (includes migration)
		if(spp_lookup$species_pooled[s] %in% c("Chum", "Pink") == FALSE) { # No rearing distributions for pink and chum
		  
		  # EXCEPTION: Harrison upstream-migrating late. This CU boundary includes spawning habitat but not the rearing lake. Use Harrison downstream migrating CU boundary to determine rearing habitat.
		  if(cu_list$cu_name_pse[which(cu_list$species_pooled == spp_lookup$species_pooled[s])[i]] == "Harrison-Upstream Migrating-Late"){
		    
		    dum.rear <- st_intersects(rearing.s, cu_boundary[cu_boundary$cuname == "Harrison-Downstream Migrating-Late", ], sparse = FALSE)
		    
		  } else { # If any other CU, use CU boubndary
		    
		    dum.rear <- st_intersects(rearing.s, cu_boundary.i, sparse = FALSE)
		  }
		  
		  rearing.i <- rearing.s[which(dum.rear == TRUE),]
		  
		  intrscts_rear <- st_intersects(rearing.i, grid_polys, sparse = FALSE)
		  incl_rear <- which(apply(intrscts_rear, 2, sum) > 0)
		}
		
		# Assign grid cells to each life stage
		I <- I + 1
		PCIC_incl[[I, "adult_migration"]] <- incl_mig
		PCIC_incl[[I, "spawning"]] <- incl_spawn
		PCIC_incl[[I, "eggs_alevin"]] <- incl_spawn
		if(spp_lookup$species_pooled[s] %in% c("Chum", "Pink") == FALSE) {
		  PCIC_incl[[I, "fw_rearing"]] <- incl_rear
		}
		
		if(I == 1){
			PCIC_incl_overall <- unique(c(incl_mig, incl_spawn, incl_rear))
		} else {
			PCIC_incl_overall <- unique(c(PCIC_incl_overall, incl_mig, incl_spawn, incl_rear))
		}
		
		# Plot
		pdf(file = paste0("freshwater/output/freshwater_distribution_fraser_", spp_lookup[s,1], "_", cuid.s[i], "_", Sys.Date() ".pdf"), width = 8, height = 10)
		# png(paste0("freshwater/output/freshwater_distribution_fraser_", spp_lookup[s,1], "_", cuid.s[i], ".pdf"))
		
		plot(st_geometry(grid_polys[unique(c(incl_mig, incl_spawn, incl_rear)),]), col = NA, border = NA, axes = TRUE) 
		plot(st_geometry(cu_boundary.i), col = grey(0.8), border = NA, add = TRUE)
		
		plot(st_geometry(grid_polys[incl_mig,]), col = paste0(cols[3], 50), border = cols[3], lwd = 0.5, add = TRUE) 
		plot(st_geometry(grid_polys[incl_spawn,]), col = paste0(cols[2], 50), border = cols[2], lwd = 0.5, add = TRUE)
		
		# Rearing (not for pink and chum)
		if(spp_lookup$species_pooled[s] %in% c("Chum", "Pink") == FALSE) {
		  plot(st_geometry(grid_polys[incl_rear,]), col = paste0(cols[1], 50), border = cols[1], lwd = 0.5, add = TRUE)
		  plot(st_geometry(rearing.i), col = cols[1], add = TRUE, lwd = 3)
		} 
		
		plot(st_geometry(shoreline), add = TRUE)
		plot(st_geometry(spawning.i), col = cols[2], add = TRUE, lwd = 2)
		plot(st_geometry(mig_paths.i), col = cols[3], add = TRUE, lwd = 0.8)
		
		legend("topright", col = c(grey(0.8), cols), lwd = c(10, 3, 2, 0.8), c("CU boundary", "Rearing", "Spawning", "Migration"), bty = "n")
		mtext(side = 3, adj = 0, paste0(cu_list$cu_name_pse[cu_list$cuid == cuid.s[i]], " ", cu_list$species_pooled[cu_list$cuid == cuid.s[i]]), font = 2)
		dev.off()
		
		print(".")
	} # end CU i
	print(paste0("End ", spp_lookup$species_pooled[s]))

	} # end species s



saveRDS(PCIC_incl, file = "freshwater/output/PCIC_incl.rds")
write.csv(PCIC_incl_overall, file = "freshwater/output/PCIC_incl_overall.csv", row.names = FALSE)
