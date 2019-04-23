
png(filename = "RRfigure.png",width = 5, height = 5, units = "cm", res = 150)
par(mar = c(0,0,0,0))
plot(0, type = "n",
     xlim = c(0,1), ylim = c(0,1),
     ylab = "", xlab = "",
     xaxt = "n", yaxt = "n", bty = "n")
rect(0,0,1,1, col = "grey80", border = NA)
rect(0,0,.75,.75, col = "grey60", border = NA)
rect(0,0,.75*.4,.75*.4, col = "grey40", border = NA)
text(.0,.85,"Target\npopulation", pos = 4, font = 2, cex = .75,offset = 0.25)
text(.0,.5,"Source\npopulation", pos = 4, font = 2, cex = .75,offset = 0.25)
text(.0,.15,"Study\nsample", pos = 4, col = "white", font = 2, cex = .75,offset = 0.25)
rect(0,0,1,1)

dev.off()