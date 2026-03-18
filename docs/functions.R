plot_cu <- function(
   CU = "Chinook: Upper Fraser River (Spring 5-2)"
    ){
  
  # Select cU
  i <- which(species_cu_name == CU)
  
  par(mfrow = c(2, 4), mar = c(3, 3, 1, 0), oma = c(1,3,4,1))

for(k in 1:2){ # For temperature and low flow
  for(j in c(3, 4, 1, 2)){
    
    plot(1:4, fw_output_summary[i, j, k, , 1, 1], "n", ylim = range(fw_output_summary[, j, k, , , 1]), xlab = "", xaxt = "n", las = 1, bty = "l", ylab = "")
    axis(side = 1, at = c(1:4), c("hist", "early", "mid", "late"))
    for(ii in 1:n.CUs){
      lines(1:4, fw_output_summary[ii, j, k, , 1, 1], col = grey(0.8))
      lines(1:4, fw_output_summary[ii, j, k, , 2, 1], col = grey(0.8))
    }
    
    for(r in 1:2){ # RCP4.5 or 8.5
      polygon(
        x = c(1:4, 4:1), 
        y = c(fw_output_summary[i, j, k, , r, 2], rev(fw_output_summary[i, j, k, , r, 3])),
        col = paste0(cols[c(3,1)[r]], 50),
        border = NA)
      points(1:4, fw_output_summary[i, j, k, , r, 1], "o", pch = 19, lwd = 2, col = cols[c(3,1)[r]])
    } # end rcp
    
    if(k == 1){
      mtext(side = 3, adj = 0, line = 1, paste0(c("c", "d", "a", "b")[j], ") ", c("Adult migration", "Spawning", "Incubation", "Freshwater rearing")[j]), cex = 0.8)
    }
    
    if(j == 3){
      mtext(side = 2, paste0("Days ", c("above\ntemp.", "below\nflow")[k], " threshold"), line = 3, , cex = 0.8)
    }
    
  } # end stage
  
  # mtext(side = 3, adj = 0, outer = TRUE, species_cu_name[i], line = 1)
  
  if(k == 1){# Add legend
    u <- par('usr')
    legend(-8, u[4] + (u[4] - u[3])*0.4, xpd = NA, lwd = 1, col = c(grey(0.8)), legend = c("All CUs, all scenarios"), bty = "n")
    legend(-3, u[4] + (u[4] - u[3])*0.4, xpd = NA, ncol = 2, lwd = 2, pch = 19, col = cols[c(3,1)], legend = c("RCP 4.5", "RCP 8.5"), bty = "n")
  }
} # end k
mtext(side = 1, outer = TRUE, "Time period", line = -0.5, cex = 0.8)

} # end function