# Code to compare optimal temperature ranges (topt) for migrating sockeye salmon
library(dplyr)
topt1 <- read.csv("1-exposure/data/optimum-temps_BC2001.csv") %>% subset(Species == "Sockeye" & Life.stage == "Migration")

topt2 <- read.csv("1-exposure/data/topt_Eliason.csv")


popCol <- c("#384E9D", "#EF1F25", "#EC7E21", "#7C2781", "#0B803F", "#6ACADC", "#FF40FF")

quartz()
par(mar = c(4,10,2,2))
plot(c(7, 26), c(1,7), "n", xlab = "Temperature (˚C)", ylab = NA, yaxt = "n", bty = "n")
axis(side = 2, at = 7:1, labels = topt2$population, las = 1)
polygon(x = c(topt1$min, topt1$min, topt1$max, topt1$max), y = c(0, 8, 8, 0), col = grey(0.8), border = NA)
segments(x0 = topt2$min, x1 = topt2$max, y0 = 7:1, y1 = 7:1, col = popCol, lwd = 10)
