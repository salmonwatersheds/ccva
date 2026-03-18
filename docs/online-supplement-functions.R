#------------------------------------------------------------------------------
# Conservation Units
#------------------------------------------------------------------------------

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

# Create short name for labelling
cu_list$cu_name_short <- cu_list$cu_name_pse
cu_list$cu_name_short <- gsub("River", "", cu_list$cu_name_short)
cu_list$cu_name_short <- gsub("Middle", "Mid", cu_list$cu_name_short)
cu_list$cu_name_short <- gsub("Lower", "Low", cu_list$cu_name_short)
cu_list$cu_name_short <- gsub("Upper", "Up", cu_list$cu_name_short)
cu_list$cu_name_short <- gsub("  ", " ", cu_list$cu_name_short)
cu_list$cu_name_short <- gsub(" -", "-", cu_list$cu_name_short)
cu_list$cu_name_short <- gsub("Migrating", "Mig", cu_list$cu_name_short)
cu_list$cu_name_short <- gsub("(cyclic)", "cycl", cu_list$cu_name_short)
cu_list$cu_name_short <- gsub("Early Summer", "ES", cu_list$cu_name_short)
cu_list$cu_name_short <- gsub("South", "S", cu_list$cu_name_short)
cu_list$cu_name_short <- gsub("North", "N", cu_list$cu_name_short)
cu_list$cu_name_short <- gsub("Creek", "Cr", cu_list$cu_name_short)

#------------------------------------------------------------------------------
# Timing
#------------------------------------------------------------------------------

# Use timing data that shortens migration window for adult lake-type sockeye to the time taken to migrate to the rearing lake, assuming that salmon holding in the lake can moderate temperatures by moving to deeper, cooler water in the lake.
timing <- read.csv("data/timing/timing-fraser.csv")
timing_days <- read.csv("data/timing/timing-fraser_days.csv")

#------------------------------------------------------------------------------
# Stages
#------------------------------------------------------------------------------

stages <- c("eggs_alevin", "fw_rearing", "early_marine", "marine_rearing", "adult_migration", "spawning") 

#------------------------------------------------------------------------------
# List of GCMs being considered
#------------------------------------------------------------------------------
gcms <- c("ACCESS1-0", "CanESM2", "CCSM4", "CNRM-CM5", "HadGEM2-ES", "MPI-ESM-LR")

#------------------------------------------------------------------------------
# Load output datasets
#------------------------------------------------------------------------------

out.dat_all <- read.csv(paste0("output/ignore/exposure-raw.csv")) # Raw among models
exposure_stage <- read.csv("output/ignore/exposure_stages.csv") # Summarised across variables for each stage
exposure <- read.csv("output/exposure_overall.csv") # Summarised across stages

#------------------------------------------------------------------------------
# Order of CUs for stage exposure figures
#------------------------------------------------------------------------------

cu_order <- c(
  # Chinook (Summer 4-1)
  313, 315, 307,
  # Chinook (Fall 4-1)
  303, 302,
 # Chinook (Spring 4-2)
 317,
   # Chinook (Spring 5-2)
  312, 318, 308, 310, 304,
 # Chinook (Summer 4-2)
 316,
 # Chinook (Summer 5-2)
  319, 314, 311, 305, 306,
  # Chinook (Fall 5-2)
  309,
  # Chinook (no brackets)
  333, 
 # Chum
 701,
 # Coho - interior to coast
 749, 709, 708, 707, 704, 705, 750, 906,
 # Pink
 710,
 # Sockeye - Early Stuart
 732,
 # SEL - Early Summer
 735, 727, 740, 752, 751, 738, 720, 734, 719,718, 711, 715,
 # SEL - Summer
 725, 721, 728, 731,
# SEL - Late
 739, 716, 714, 713, 729, 712,
 # SER
 745,742,
 # Steelhead (Summer)
 780, 781, 783, 784, 
# Steelhead (winter)
782, 785
 
)

###############################################################################
# Figure functions
###############################################################################

species_pch <- c(Chinook = 22, Chum = 24, Coho = 21, Pink = 9, Sockeye = 23, Steelhead = 25)

#------------------------------------------------------------------------------
# Overall exposure
#------------------------------------------------------------------------------
plot.overallExposure <- function(emissions = "rcp45", period = "mid", s = 1, order.type = "lifeHistory", cu.num = FALSE){
  # quartz(width = 7.2, height = 8.5, pointsize = 10)
  r <- match(emissions, c("rcp45", "rcp85"))
 p <- match(period, c("hist", "early", "mid", "late"))
 
 rcp <- c("rcp45", "rcp85")[r]
 rcp_display <- c("RCP 4.5", "RCP 8.5")[r]
 period_display <- c("Historical (1970-1999)", "Early century (2010-2039)", "Mid century (2040-2069)", "Late century (2070-2099)")[p]
 
 ind <- which(exposure$rcp == rcp & exposure$period == c("hist", "early", "mid", "late")[p])
 var <- which(names(exposure) == c("exp_simple", "exp_weighted")[s])
 exp_summary0 <- cbind(
   tapply(exposure[ind, var], exposure$cuid[ind], median),
   tapply(exposure[ind, var], exposure$cuid[ind], min),
   tapply(exposure[ind, var], exposure$cuid[ind], max)
 )
 
 
 if(order.type == "lifeHistory"){ # Order CUs by species and life-history type
   exp_summary <- exp_summary0[match(cu_order, as.numeric(row.names(exp_summary0))), ]
   cu_list.plot <- cu_list[match(cu_order, cu_list$cuid),]
 } else if(order.type == "decreasingExposure"){ # Alternative: order by decreasing median exposure
   cu_order.exp <- as.numeric(rownames(exp_summary0))[order(exp_summary0[, 1], decreasing = TRUE)]
   exp_summary <- exp_summary0[match(cu_order.exp, as.numeric(row.names(exp_summary0))), ]
   cu_list.plot <- cu_list[match(cu_order.exp, cu_list$cuid),]
 }
   
 # Create palette for overall exposure, over range of median for rcp and period
 col_palette <- grey.colors(n = n+1, start = 0.99, end = 0.01)
 col_levels_overall <- seq(min(exp_summary[,1]), max(exp_summary[,1]), length.out = n)
 
 par(mar = c(4,16,2,1))
 plot(exp_summary[, 1], c(60:1), "n", xlim = range(exp_summary), xlab = "", yaxt = "n", ylab = "", bty = "l", yaxs = "i", ylim = c(0.5, 60.5))
 u <- par('usr')
 
 mtext(side = 1, line = 2, "Overall exposure")
 
 # Add horizontal and vertical guides
 abline(h = seq(1, 62, 2) - 0.5, lty = 3, col = grey(0.8), xpd = NA)
 abline(v = seq(0, 1, 0.05), lty = 3, col = grey(0.8))
 
 # Axis
 axis(side = 1, at = seq(0, 1, 0.05), tck = -0.01, labels = FALSE)
 
text(u[1], 61.5, paste("Conservation Unit", " "), cex = 0.8,adj = 1, xpd = NA)
 text(u[1], c(60:1), paste(cu_list.plot$cu_name_short, " "), cex = 0.8, col = species_cols[cu_list.plot$species_pooled], adj = 1, xpd = NA)
 
 segments(x0 = exp_summary[, 2], x1 = exp_summary[, 3], y0 = c(60:1), y1 = c(60:1), col = species_cols[cu_list.plot$species_pooled], lwd = 2)
 
 # legend(par('usr')[1] + 0.05*(par('usr')[2] - par('usr')[1]), 60, col = species_cols, pch = 19, pt.cex = 1.5, lwd = 2, names(species_cols), bg = "white")
 # legend("bottomright", col = species_cols, pch = 19, pt.cex = 1.5, lwd = 2, names(species_cols), bty = "n")
 
 # Layer on points

 if(cu.num == TRUE){
   points(exp_summary[, 1][which(row.names(exp_summary) == "710")], c(60:1)[which(row.names(exp_summary) == "710")], bg = col_palette[findInterval(exp_summary[, 1], col_levels_overall)][which(row.names(exp_summary) == "710")], lwd = NA, col= NA, cex = 2.5, pch = 23) # Pink background
   points(exp_summary[, 1], c(60:1), col = species_cols[cu_list.plot$species_pooled], pch = species_pch[cu_list.plot$species_pooled], cex = 2.5, bg = col_palette[findInterval(exp_summary[, 1], col_levels_overall)], lwd = 2, xpd = NA)
   text(exp_summary[, 1], c(60:1), c(1:60), cex = 0.7, font = 2)
   text(exp_summary[31, 1], 30, 31, cex = 0.7, font = 2, col = "white")
   
 } else {
   points(exp_summary[, 1][which(row.names(exp_summary) == "710")], c(60:1)[which(row.names(exp_summary) == "710")], bg = col_palette[findInterval(exp_summary[, 1], col_levels_overall)][which(row.names(exp_summary) == "710")], lwd = NA, col= NA, cex = 1.5, pch = 23) # Pink background
   points(exp_summary[, 1], c(60:1), col = species_cols[cu_list.plot$species_pooled], pch = species_pch[cu_list.plot$species_pooled], cex = 1.5, bg = col_palette[findInterval(exp_summary[, 1], col_levels_overall)], lwd = 2, xpd = NA)
 }

 mtext(side = 3, line = -1, paste(period_display, "under", rcp_display, sep = " "), outer = TRUE)
 
 if(order.type == "lifeHistory"){
   sp_start <- c(1, 20, 21, 29, 30, 55, 700)
   xx <- u[1] - 0.55*(u[2] - u[1])
   xxx <- u[1] - c(0.56, 0.57)[c(1,2,1,2,1,1)]*(u[2] - u[1])
   for(s in 1:6){
     segments(x0 = xx, x1 = xx,
              y0 = 60 - sp_start[s] + 1.2,
              y1 = max(0.8, 60 - sp_start[s+1] + 1.7),
              col = species_cols[s],
              xpd = NA, lwd = 2)
     text(xxx[s], mean(c(60 - sp_start[s] + 1.2, max(0.8, 60 - sp_start[s+1] + 1.7))) - 1, c("Chinook", "Chum", "Coho", "Pink", "Sockeye", "Steelhead")[s], srt = c(90, 0, 90, 0, 90, 90)[s], xpd = NA, font = 2, col = species_cols[s], pos = 3)
   }
 }
  
  }

#------------------------------------------------------------------------------
# Life stage exposure
#------------------------------------------------------------------------------
# Look at mean and range across models
stage_summary <- array(NA, 
                       dim = c(n.CUs, length(stages), 2, 2, 4), 
                       dimnames = list(cuid, stages, NULL, c("rcp45", "rcp85"), c("hist", "early", "mid", "late")))

for(i in 1:n.CUs){
  for(j in 1:length(stages)){
    for(r in 1:2){
      for(p in 1:4){
        for(v in 1:2){
      ind <- which(out.dat_all$cuid == cuid[i] & out.dat_all$stage == stages[j] & out.dat_all$rcp == c("rcp45", "rcp85")[r] & out.dat_all$period == c("hist", "early", "mid", "late")[p] & out.dat_all$variable == list(c("stream_temp", "stream_flow"), c("SST", "SSS"))[[c(1,1,2,2,1,1)[j]]][v])
      
      # Simple
      stage_summary[i, j, v, r, p] <- median(out.dat_all$median[ind])
      
      }
    }
  }
  }
} # end cu i

plot.stageExposure <- function(emissions = "rcp45", period = "mid"){
  r <- match(emissions, c("rcp45", "rcp85"))
  p <- match(period, c("hist", "early", "mid", "late"))

  rcp <- c("rcp45", "rcp85")[r]
  rcp_display <- c("RCP 4.5", "RCP 8.5")[r]
  period_display <- c("Historical (1970-1999)", "Early century (2010-2039)", "Mid century (2040-2069)", "Late century (2070-2099)")[p]
  
  # order of species
  #----
  # Create o order to match order above? or just set o <- 1:60 for default ordering of CUs by species
  ind <- which(exposure$rcp == "rcp45" & exposure$period == "mid")
  var <- which(names(exposure) == "exp_simple")
  exp_summary0 <- cbind(
    tapply(exposure[ind, var], exposure$cuid[ind], median),
    tapply(exposure[ind, var], exposure$cuid[ind], min),
    tapply(exposure[ind, var], exposure$cuid[ind], max)
  )

  # re-order by cu_list cuids
  dum.cuid <- as.numeric(rownames(exp_summary0))
  exp_summary <- exp_summary0[match(cuid, dum.cuid), ]

  # Check:
  if(sum(as.numeric(rownames(exp_summary)) - cuid) > 0){  # Should be zero
    stop("CUID order doesn't match!!")
  }
  # o <- order(exp_summary[, 1], decreasing = TRUE)
  # o <- c(
  #   # Chinook
  #   order(exp_summary[which(cu_list$species_pooled == "Chinook"), 1], decreasing = TRUE),
  #   # Chum
  #   length(which(cu_list$species_pooled %in% c("Chinook", "Chum"))),
  #   # Coho
  #   length(which(cu_list$species_pooled %in% c("Chinook", "Chum"))) + order(exp_summary[which(cu_list$species_pooled == "Coho"), 1], decreasing = TRUE),
  #   # pink     
  #   length(which(cu_list$species_pooled %in% c("Chinook", "Chum", "Coho", "Pink"))),
  #   # Sockeye     
  #   length(which(cu_list$species_pooled %in% c("Chinook", "Chum", "Coho", "Pink"))) + order(exp_summary[which(cu_list$species_pooled == "Sockeye"), 1], decreasing = TRUE),
  #   # Steelhead     
  #   length(which(cu_list$species_pooled %in% c("Chinook", "Chum", "Coho", "Pink", "Sockeye"))) + order(exp_summary[which(cu_list$species_pooled == "Steelhead"), 1], decreasing = TRUE)
  # ) 
  
  # o <- c(
  #   order(apply(stage_summary[which(cu_list$species_pooled == "Chinook"), , , r, p], 1, mean), decreasing = TRUE),
  #   length(which(cu_list$species_pooled %in% c("Chinook", "Chum"))),
  #   length(which(cu_list$species_pooled %in% c("Chinook", "Chum"))) + order(apply(stage_summary[which(cu_list$species_pooled == "Coho"), , , r, p], 1, mean), decreasing = TRUE),
  #   length(which(cu_list$species_pooled %in% c("Chinook", "Chum", "Coho", "Pink"))),
  #   length(which(cu_list$species_pooled %in% c("Chinook", "Chum", "Coho", "Pink"))) + order(apply(stage_summary[which(cu_list$species_pooled == "Sockeye"), , , r, p], 1, mean), decreasing = TRUE),
  #   length(which(cu_list$species_pooled %in% c("Chinook", "Chum", "Coho", "Pink", "Sockeye"))) + order(apply(stage_summary[which(cu_list$species_pooled == "Steelhead"), , , r, p], 1, mean), decreasing = TRUE)
  # )
    
  o <- match(cu_order, as.numeric(dimnames(stage_summary)[[1]]))
    
  
  # quartz(width = 7.1, height = 9.5, pointsize = 10)
  par(mar = c(3.5,16,7.5,1))
  plot(c(0,6), c(60,1), "n", xlab = "", yaxt = "n", xaxs = "i", ylab = "", bty = "l", yaxs = "i", ylim = c(0.5, 60.5), xaxt = "n")
  mtext(side = 1, line = 2.5, "Life stages")
  # Layer on points
  for(j in 1:6){
    z <- stage_summary[o, j, , r, p]
    for(i in 60:1){
      for(v in 1:2){
        polygon(x = (j - 1) + 0.5*(v-1) + c(0, 0, 0.5, 0.5),
                y = (i - 0.5) + c(0, 1, 1, 0),
                border = col_palette[findInterval(z[60 - i + 1, v], col_levels)],
                col = col_palette[findInterval(z[60 - i + 1, v], col_levels)]
        )
      }
      
    } # end i
  } # end j
  # Add horizontal and vertical guides
  abline(h = seq(1, 62, 2) - 0.5, lty = 3, col = grey(0.8), xpd = NA)
  abline(v = seq(0.5, 6.5, 1), lty = 3, col = grey(0.8))
  abline(v = seq(1, 6, 1), lwd = 2)
  
  
  # axis(side = 2, at = c(60:1), labels = FALSE)
  text(par('usr')[1], 61.5, paste("Conservation Unit", " "), cex = 0.9, adj = 1, xpd = NA)
  text(par('usr')[1], c(60:1), paste(cu_list$cu_name_short[o], " "), cex = 0.8, col = species_cols[cu_list$species_pooled[o]], adj = 1, xpd = NA)
  
  # Axis
  axis(side = 1, at = seq(0, 6, 1), labels = FALSE)
  text(seq(0.5, 5.5, 1), -1, c("incubation", "freshwater\nrearing", "early\nmarine", "marine\nrearing", "adult\nmigration", "spawning"), xpd = NA, cex = 0.9)
  
  # axis(side = 3, at = seq(0.5, 5.5, 1), labels = FALSE)
  # text(seq(0.5, 5.5, 1), 63, c("Incubation", "FW\nrear.", "Early\nmarine", "Marine\nrearing", "Adult\nFW mig.", "Spawn"), xpd = NA)
  axis(side = 3, at = seq(0, 6, 0.5), labels = FALSE)
  # text(seq(0.25, 5.75, 0.5), 62.5, c(rep(c("FWT", "FWF"), 2), rep(c("SST", "SSS"), 2), rep(c("FWT", "FWF"), 2)), srt = 90, xpd = NA)
  text(seq(0.25, 5.75, 0.5), 62.2, c(rep(c("temp", "flow"), 2), rep(c("SST", "sal"), 2), rep(c("temp", "flow"), 2)), srt = 90, xpd = NA, cex = 0.9)
  
  text(par('usr')[1], c(60:1), paste(cu_list$cu_name_short[o], " "), cex = 0.8, col = species_cols[cu_list$species_pooled[o]], adj = 1, xpd = NA)
  
  # Legend
  xw <- c(0, 6)
  yw <- c(66.5, 68)
  polygon(x = xw[c(1,2,2,1)], y = yw[c(1,1,2,2)], xpd= NA, col = "#FFFFFF")
  for(k in 1:n){
    polygon(x = xw[1] + (k-1)*diff(xw)/n + c(0, diff(xw)/n, diff(xw)/n, 0),
            y = yw[c(1,1,2,2)],
            col = col_palette[k],
            border = NA,
            xpd = NA)
  }
  segments(x0 = seq(xw[1], xw[2], length.out = 5), x1 = seq(xw[1], xw[2], length.out = 5), y0 = yw[1]-0.5, y1 = yw[2], xpd = NA)
  
  text(seq(xw[1], xw[2], length.out = 5), yw[1]-1, seq(0, 1, length.out = 5), xpd = NA)
  text(xw[1], yw[1], pos = 2, "Proportion of time above\nor below threshold", xpd = NA)
  mtext(side = 3, line = -1, paste(period_display, "under", rcp_display, sep = " "), outer = TRUE)
  
  sp_start <- c(1, 20, 21, 29, 30, 55, 700)
  u <- par('usr')
  xx <- u[1] - 0.55*(u[2] - u[1])
  xxx <- u[1] - c(0.56, 0.57)[c(1,2,1,2,1,1)]*(u[2] - u[1])
  for(s in 1:6){
    segments(x0 = xx, x1 = xx,
             y0 = 60 - sp_start[s] + 1.2,
             y1 = max(0.8, 60 - sp_start[s+1] + 1.7),
             col = species_cols[s],
             xpd = NA, lwd = 2)
    text(xxx[s], mean(c(60 - sp_start[s] + 1.2, max(0.8, 60 - sp_start[s+1] + 1.7))) - 1, c("Chinook", "Chum", "Coho", "Pink", "Sockeye", "Steelhead")[s], srt = c(90, 0, 90, 0, 90, 90)[s], xpd = NA, font = 2, col = species_cols[s], pos = 3)
  }
  
}

#------------------------------------------------------------------------------
# Last figure: Change over periods
#------------------------------------------------------------------------------

plot.overPeriods <- function(rcp_scen = "rcp45"){
  
 stage.names <- c("Incubation", "Freshwater rearing", "Early marine", "Marine rearing", "Adult migration", 'Spawning')
  var.names <- cbind(c("stream temperature", "stream flow"), c("SST", "salinity"))
  
  par(mfrow = c(3, 4), oma = c(1,10,1,0))
  for(j in 1:6){
    for(v in 1:2){
      if (v == 1){
        par(mar = c(4, 1, 3, 0))
      } else {
        par(mar = c(4, 0, 3, 1))
      }
      
      var <- list(c("stream_temp", "stream_flow")[v], c("SST", "SSS")[v])[[c(1, 1, 2, 2, 1, 1)[j]]]
      plot(c(1,4), c(1,7), "n", yaxs = "i", xaxs = "i", xlim = c(0.8, 4.2), ylim = c(1, 7), bty = "n", xaxt = "n", xlab = "", ylab = "", yaxt = "n")
      axis(side = 1, at = c(1:4), c("H", "E", "M", "L"))
      abline(h = 1:6, lty = 3, col = grey(0.8))
      abline(v = 1:4, lty = 3, col = grey(0.8))
      abline(v = 3, lty = 3)
      mtext(side = 3, var.names[v, c(1,1,2,2,1,1)[j]], line = 0.2)
      if(v == 1) text(4.2, 7.9, xpd = NA, stage.names[j], cex = 1.6)
      if (j %in% c(1,3,5) & v == 1){
        mtext(side = 2, line = 8, "Exposure")
        for(s in 1:6){
          axis(side = 2, at = 7 - c(s, s - 1), label = c(0,1), las = 1, line = rep(c(-0.5, 1.5), 3)[s], col = species_cols[s])
          text(rep(c(0.6, 0.1), 3)[s], 7 - s + 0.5, names(species_cols)[s], col = species_cols[s], xpd = NA, adj = 1, cex = 1.3)
        }}
      for(s in 1:6){
        spql <- list("CK", "CM", "CO" ,"PKO", c("SEL", "SER"), "SH")[[s]]
        out.dat.js <- out.dat_all %>% 
          filter(rcp == rcp_scen, species %in% spql, stage == stages[j], variable == var)
        for(i in unique(out.dat.js$cuid)){
          out.dat.jsi <- out.dat.js %>% filter(cuid == i)
          y <- cbind(median = tapply(out.dat.jsi$median, out.dat.jsi$period, median),
                     min = tapply(out.dat.jsi$median, out.dat.jsi$period, min),
                     max = tapply(out.dat.jsi$median, out.dat.jsi$period, max))
          
          lines(1:4, 7 - s + y[,1], col = species_cols[s])
          polygon(x = c(1:4, 4:1),
                  y = c(7 - s + y[,2], rev(7 - s + y[,3])),
                  col = paste0(species_cols[s], 20), border = NA)
        } # end cu i 
      }# end s
      
    } # end v
  } # end j
  
  mtext(side = 1, line = 0, outer = TRUE, "Period")
}