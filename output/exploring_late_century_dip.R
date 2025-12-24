cols <- wes_palette("Darjeeling1")
cols_period <- cols[c(5, 3, 4, 1)]
names(cols_period) <- c("hist", "early", "mid", "late")

# Salinity
source("marine/code/marine-functions.R")

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

quartz(width = 7, height = 8, pointsize = 10)
par(mfrow = c(3,2), mar = c(3, 3, 2, 1), oma = c(2,2,0,0))
for(m in 1:6){
rcp <- 85 
SSS <- loadNOAA(
  variable = "SSS", # Which variable to load? One of "PP", "pH", "SSS", "SST"
  model = gcms$modelName[m], # Which GCM? 
  rcp = rcp, # Which emissions scenario? Options 45 or 85
  spat_id = spat_EM 
)


# Check dimnames: this is how output will be organized
grid.ref <- as.numeric(dimnames(SSS)[[1]])

dim(SSS)

# Plot

plot(months, SSS[1,], "n", xlim = c(2.5,9.5), xaxs = "i")
abline(v = seq(3.5, 8.5, 1))
for(p in 1:4){
  P <- c("hist", "early", "mid", "late")[p]
  points(months[which(period == P)]+(p-2.5)/4, SSS[1,which(period == P)], col = cols_period[P])
  lines(1:12, tapply(SSS[1,which(period == P)], months[which(period == P)], median), lwd = 2, col = cols_period[P])
}
abline(h = quantile(SSS[1,which(period == "hist")], 0.025), lty= 2)
mtext(side = 3, gcms$modelName[m])
} # end m

mtext(side = 1, outer = TRUE, "Month")
mtext(side = 2, outer = TRUE, "Sea surface salinity (ppt)")
