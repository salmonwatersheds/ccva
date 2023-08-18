###############################################################################
# Functions to load and summarize data relevant to exposure indicators
# 
# Contact: Stephanie Peacock <speacock@psf.ca>
# Date: April 17, 2023
###############################################################################


#------------------------------------------------------------------------------
# Function to calculate 7-day average 
#------------------------------------------------------------------------------

calc_7DADM <- function(
    Tmax, # time-series of maximum daily temperatures
    date # associated dates
){
  
  if(sum(unique(diff(date)) != 1) > 0){
    stop("Not continuous time series")
    
  } else {
    
    n <- length(Tmax)
    
    dum <- rbind(
      c(rep(NA, 6), Tmax[1]),
      c(rep(NA, 5), Tmax[1:2]),
      c(rep(NA, 4), Tmax[1:3]),
      c(rep(NA, 3), Tmax[1:4]),
      c(rep(NA, 2), Tmax[1:5]),
      c(NA, Tmax[1:6]),
      cbind(
        Tmax[c(1:(n - 6))],
        Tmax[c(2:(n - 6 + 1))],
        Tmax[c(3:(n - 6 + 2))],
        Tmax[c(4:(n - 6 + 3))],
        Tmax[c(5:(n - 6 + 4))],
        Tmax[c(6:(n - 6 + 5))],
        Tmax[c(7:n)])
    )
    
    T7 <- apply(dum, 1, mean, na.rm = TRUE)
    return(T7)
  }
}

#------------------------------------------------------------------------------
# Function to calculate rolling 7-day maximum
#------------------------------------------------------------------------------
weeklyStat <- function(x, stat){
  n <- length(x)
  X <- rbind(
    x, 
    c(NA, x[1:(n-1)]),
    c(rep(NA, 2), x[c(1:(n - 2))]),
    c(rep(NA, 3), x[c(1:(n - 3))]),
    c(rep(NA, 4), x[c(1:(n - 4))]),
    c(rep(NA, 5), x[c(1:(n - 5))]),
    c(rep(NA, 6), x[c(1:(n - 6))])
  )    
  if(stat == "max"){
    xx <- apply(X, 2, max, na.rm = TRUE)
  } else if(stat == "mean"){
    xx <- apply(X, 2, mean, na.rm = TRUE)
  }

  return(xx)
    
}
