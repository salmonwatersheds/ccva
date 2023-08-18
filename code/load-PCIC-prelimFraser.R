###############################################################################
# Script to load PCIC Preliminary model output
# 
# Contact: Stephanie Peacock <speacock@psf.ca>
# Date: April 17, 2023
###############################################################################

link_dat <- "data/raw-data/PCIC/hydro_model_out/Fraser-prelim/"

#-----------------------------------------------------------------------------
# Open
#-----------------------------------------------------------------------------

z <- nc_open(paste0(link_dat, "MOUTH_CanESM2_r1i1p1_rcp85_sel-001.nc"))

# # Explore dataset. What are the variables?
# print(paste("File",z$filename,"contains",z$nvars,"variables"))
# for( i in 1:z$nvars ) {
#   v <- z$var[[i]]
#   print(paste("Here is information on variable number",i))
#   print(paste("   Name: ",v$name))
#   print(paste("   Units:",v$units))
#   print(paste("   Missing value:",v$missval))
#   print(paste("   # dimensions :",v$ndims))
#   print(paste("   Variable size:",v$varsize))
# }

# Extract temperature (Ts) or flow (Qs) variables
Ts <- ncvar_get(z, "Ts")
Qs <- ncvar_get(z, "Qs")

# Create time variable
date <- as.Date(c(0:56612), origin = "1945-01-01")

rm(z)

#-----------------------------------------------------------------------------
# Load spatial grid: x 0.0625 degrees
#-----------------------------------------------------------------------------

# Subset to remove duplicate center points with V7 == 2
spat_grid0 <- read.csv(paste0(link_dat, "MOUTH.Spat.csv"), header = FALSE, col.names = c("V1", "V2", "V3", "V4", "lat", "lon", "V7", "V8", "V9")) %>% subset(V7 == 1)

# V2 is the unique identifier
# Note that in the prelim data there were repeat lat/lon points where two rivers meet
# Choose the one with the highest flow (on average?)
spat_grid0$avgQs <- apply(Qs, 1, mean)

spat_grid0$latlon <- paste(spat_grid0$lat, spat_grid0$lon, sep = "-")
spat_grid0$keep <- rep(0, nrow(spat_grid0))
for(i in 1:nrow(spat_grid0)){
  # Find indices of grid points with same coordinates
  ind <- which(spat_grid0$latlon == spat_grid0$latlon[i])
  if(length(ind) > 1){
    avgQs.i <- spat_grid0$avgQ[ind]
    if(spat_grid0$avgQ[i] == max(avgQs.i)){
      # If that point has the max flow of all overlapping points, then keep it
      spat_grid0$keep[i] <- 1
    }
  } else {
    # If there is only one of those points, keep it.
    spat_grid0$keep[i] <- 1
  }
}
# sum(spat_grid0$keep) # 8311 grid cells
spat_grid <- subset(spat_grid0, spat_grid0$keep == 1)

# Make grid center points

z <- data.frame(lon = as.numeric(spat_grid$lon), lat = as.numeric(spat_grid$lat))
# grid_points <- st_multipoint(z) %>% st_sfc(crs = 4269)
grid_points <- st_as_sf(z, coords = c("lon", "lat"), crs = 4269)

# Make grid polygon
n <- dim(spat_grid)[1]
y <- list(); length(y) <- n
d <- 1/16

poly_dat <- st_as_sf(data.frame(
  V2 = rep(spat_grid$V2, each = 5),
  rep(c("SW0", "NW", "NE", "SE", "SW1"), n),
  lon = c(rep(spat_grid$lon, each = 5) + rep(c(- d/2, -d/2, d/2,  d/2, -d/2), n)),
  lat = c(rep(spat_grid$lat, each = 5) + rep(c(- d/2,  d/2, d/2, -d/2, -d/2), n))
), coords = c("lon", "lat"), crs = 4269)

poly_dat %>% 
  group_by(V2) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") -> grid_polys

rm(poly_dat)

#-----------------------------------------------------------------------------
# Subset temperature and flow to kept grid cells only
#-----------------------------------------------------------------------------

Qs <- Qs[spat_grid$V2, ]
Ts <- Ts[spat_grid$V2, ]
