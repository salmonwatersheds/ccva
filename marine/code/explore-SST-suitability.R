maropt <- read.csv("data/Langan2024_SSTsuitability.csv")

sp_cols <- c(wesanderson::wes_palette("Darjeeling1"), grey(0.7))
names(sp_cols) <- c("Chinook", "Chum", "Coho", "Pink", "Sockeye", "Steelhead")

# Simply plot them all
plot(maropt$SST, maropt[, c("Chinook", "Chum", "Coho", "Pink_Even", "Sockeye", "Steelhead")[s]], "n", ylim = c(0,1))
for(s in 1:6){
  lines(maropt$SST, maropt[, c("Chinook", "Chum", "Coho", "Pink_Even", "Sockeye", "Steelhead")[s]], col = sp_col[s])
  if(s == 4) lines(maropt$SST, maropt[, "Pink_Odd"], col = sp_col[s], lty = 2)
}


par(mar = c(4, 2, 5, 2))
plot(maropt$SST, maropt$Sockeye, 'n', ylim = c(0, 6), yaxs = "i", bty = "n", xlab = "Temperature (˚C)", ylab = "", las = 1, yaxt = "n", xlim = c(0, 20), xaxs = "i")
mtext(side = 2, line = 0.5, "Thermal suitability")
abline(v = seq(0, 20, 1), col = grey(0.8), lwd = 0.5)
abline(h = c(1:5))

for(s in 1:6){
  y.s <- maropt[, c("Chinook", "Chum", "Coho", "Pink_Even", "Sockeye", "Steelhead")[s]]
  # x <- sample(maropt$SST, size = 10^6, replace = TRUE, prob = y.s)
  # HDInterval::hdi(x, credMass = 0.5)
  # lines(maropt$SST, (6 - s) + y.s, "l", col = sp_cols[s])
  polygon(x = c(maropt$SST, rev(maropt$SST)), y = c((6 - s) + y.s, rep(6 - s, length(y.s))), col = paste0(sp_cols[s], 60), border = NA)
  
  ind <- tail(which(y.s > 0.5), 1) + 1
  segments(x0 = maropt$SST[ind], y0 = 6 - s + y.s[ind], x1 = maropt$SST[ind], y1 = 6 - s, sp_cols[s])
  points(maropt$SST[ind], 6 - s + y.s[ind], pch = 21, col = sp_cols[s], bg = "white", lwd = 1.5, cex = 1.5)
  
  if(s == 4){ #Even and odd pink
    
    y.s <- maropt[, "Pink_Odd"]
    # lines(maropt$SST, (s-1) + y.s, "l", col = sp_cols[s])
    polygon(x = c(maropt$SST, rev(maropt$SST)), y = c(6 - s + y.s, rep(6 - s, length(y.s))), col = paste0(sp_cols[s], 60), border = NA, density = 30, angle = 45, lwd = 1.5)
    
    ind <- tail(which(y.s > 0.5), 1) + 1
    segments(x0 = maropt$SST[ind], y0 = 6 - s + y.s[ind], x1 = maropt$SST[ind], y1 = 6 - s, sp_cols[s], lty = 2)
    points(maropt$SST[ind], 6 - s + y.s[ind], pch = 21, col = sp_cols[s], bg = "white", lwd = 1.5, cex = 1.5, lty = 3)
    
  }
  
  points(as.numeric(fwopt$Temperature.threshold[fwopt$Species == c("Chinook", "Chum", "Coho", "Pink", "Sockeye", "Steelhead")[s]])[c(1,2,4)], rep(6 - s, 3), pch = c(4, 25, 24), lwd = 1.5, bg = "white", col = sp_cols[s], xpd = NA, cex = 1.5)
  segments(x0 = 18, x1 = 18, y0 = 6 - s + 1, y1 = 6 - s, lty = 3, lwd = 1.5, col = sp_cols[s])
}

text(0, 6 - c(1:6) + 1 - 0.2, c("Chinook", "Chum", "Coho", "Pink", "Sockeye", "Steelhead"), col = sp_cols, xpd = NA, font = 2, adj = 0)

legend(14, 7.2, pch = c(21, NA, 4, 24, 25), legend = c("Marine (core range upper)", "Adult freshwater migration", "Spawning", "Incubation", "Juvenile rearing"), pt.lwd = 1.5, xpd = NA, lty = c(NA, 3, NA, NA, NA), lwd = c(NA, 1.5, NA, NA, NA), title = "Upper temp. thresholds")

###############################################################################
topt <- read.csv("data/optimum-temps.csv")

stages <- unique(topt$stage)[c(4,2,5,3,1)]
stage_pch <- c(19, 24, 25, 22, 4)

species <- c("Chinook", "Chum", "Coho", "Pink", "Sockeye", "Steelhead")


plot(c(7, 19), c(0,6), "n", yaxs = "i", bty = "n", xlab = "Temperature (˚C)", ylab = "", las = 1, yaxt = "n", xaxs = "i")
abline(v = seq(0, 20, 1), col = grey(0.8), lwd = 0.5)
abline(h = c(1:6))

for(s in 1:6){
  ind <- which(topt$species == species[s])
  points(topt$max[ind][match(stages, topt$stage[ind])], 6 - s + rev(seq(0.15, 0.85, 0.15)), pch = stage_pch, lwd= 2, cex = 1.5, col = sp_cols[s])
  
}
text(7, 6 - c(1:6) + 1 - 0.2, c("Chinook", "Chum", "Coho", "Pink", "Sockeye", "Steelhead"), col = sp_cols, xpd = NA, font = 2, adj = 0)

legend(8, 7, pch = stage_pch, legend = c("Incubation", "FW juvenile", "Marine", "Adult FW mig.", "Spawning"), pt.lwd = 1.5, xpd = NA, ncol = 5, pt.cex = 1.5)
