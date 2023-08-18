shift <- c(0, 1:7, seq(10, 21, 3), seq(25, 40, 5))

Ts.dd_shift <- array(NA, dim = c(length(shift), 4, 31))
meanQs_shift <- array(NA, dim = c(length(shift), 4, 31))

for(j in 1:length(shift)){
  
#-----------------------------------------------------------------------------
# Calculate rolling weekly max (Ts) or mean (Qs) over spawning period
#-----------------------------------------------------------------------------

  # Calculate MAD - historical
  
  mad <- apply(Qs[incl, which(date %in% c(periods$start_num[1]:periods$end_num[1]))], 1, mean)
  
# Select dates in spawn timing, extending earlier 7 days for calculating 7-day rolling avg
periods.ind <- rbind(
  hist = which(DOY %in% c(c(timing.i$spawn_start - 7):timing.i$spawn_end) & years %in% c(1970:1999)),
  early = which(DOY %in% c(c(timing.i$spawn_start - 7):timing.i$spawn_end) & years %in% c(2010:2039)),
  mid = which(DOY %in% c(c(timing.i$spawn_start - 7):timing.i$spawn_end) & years %in% c(2040:2069)),
  late = which(DOY %in% c(c(timing.i$spawn_start - 7):timing.i$spawn_end) & years %in% c(2070:2099)))

# Create arrays
Ts.weeklyMax <- array(
  NA, 
  dim = c(4, length(incl), ncol(periods.ind) - 7), 
  dimnames = list(c("hist", "early", "mid", "late"))

for(p in 1:4){ # For each period
  Ts.weeklyMax[p, , ] <- t(apply(Ts[incl, periods.ind[p,]], 1, weeklyStat, stat = "max"))[, 8:ncol(periods.ind)]
  Qs.weeklyMean[p, , ] <- t(apply(Qs[incl, periods.ind[p,]], 1, weeklyStat, stat = "mean"))[, 8:ncol(periods.ind)]
}

# Checks in case -Inf (not getting this any more?)
Ts.weeklyMax[which(Ts.weeklyMax == -Inf)] <- NA
Qs.weeklyMean[which(Qs.weeklyMean == -Inf)] <- NA

# Average weekly max for each period 
# Calculate degree days (dd) pver thresh
thresh <- 18
Ts.weeklyMax_avg <- matrix(NA, nrow = 4, ncol = 31)
Ts.dd <- matrix(NA, nrow = 4, ncol = 31)
for(p in 1:4){
  Ts.weeklyMax_avg[p, ] <- apply(Ts.weeklyMax[p, , ], 1, mean, na.rm = TRUE)
  Ts.dd[p, ] <- apply(Ts.weeklyMax[p, , ], 1, function(x){sum(pmax(0, x - thresh))/30})
}


###############################################################################