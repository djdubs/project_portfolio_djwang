# Author: Daniel Wang
# Last update: 6/24/2026

library(spatstat)
library(sf)
library(splancs)
library(dplyr)
library(ggplot2)
library(viridis)
library(PlotTools)

wake=read.csv("projects/data/WakeCrime.csv")
shapeWake <- st_read("projects/data/Townships-shp/townships.shp")
shapeWake=st_transform(shapeWake, "EPSG:4326")

plot(shapeWake["NAME"])

# Polygon coordinates for Neuse
Neuse=shapeWake[[5]][shapeWake$NAME=="NEUSE"]
Neuse.poly=st_coordinates(Neuse)[,1:2]

# table of crimes by common occurrence
sort(table(wake$Crime_Category), decreasing = T)

# I am interested in examining the point patterns of traffic crimes on weekends vs weekdays
weekends <- c("Saturday", "Sunday")
daytime <- seq(6,17)

# respective data sets for night and day times
wake1 <- wake[wake$Reported.Day.of.Week %in% weekends &
                wake$Crime_Category=="TRAFFIC",]
wake2 <- wake[!(wake$Reported.Day.of.Week %in% weekends) &
                wake$Crime_Category=="TRAFFIC",]

# subsetting to rows in Neuse
N1 <- pip(as.points(wake1$Longitude, wake1$Latitude), Neuse.poly, out = F)
N2 <- pip(as.points(wake2$Longitude, wake2$Latitude), Neuse.poly, out = F)

nsim=5
n <- nrow(N1)
n2 <- nrow(N2)

# simple demo plot to view pt patterns
png("projects/results/sa_p1.png", width = 8, height = 5, units = "in", res=300)
par(mar = c(5.1,4.1,1,1))
plot(Neuse.poly, type="l", asp=1, cex.axis = 1, xlab=NA, ylab=NA)
title(xlab = "Longitude", ylab = "Latitude", cex.lab = 1)
points(N1,col=2, pch=1)
points(N2,col=4, pch=3)
legend("topleft",c("Weekend","Weekday"),col=c(2,4),pch=c(1,3), cex = 1)
par(mar = c(5.1,4.1,4.1,2.1))
dev.off()

# ______________________________________________________________________________
# Assessing both point patterns for spatial randomness
# ______________________________________________________________________________
# using a maximum spatial lag of 3km (0.03 degrees)

# K function for weekend point pattern

# computing k function over spatial lag
h=seq(.001,0.03,0.001)
kpts=khat(N1,Neuse.poly,h)

# k function plot
png("projects/results/sa_p2.png", width = 800, height = 600)
par(mar = c(5.1,4.1,1,1))
plot(h, kpts, type="l", lwd=1.5, ylab="K", xlab="Spatial Lag",
     main = "K-function vs. CSR curve | Weekend Group")
points(h, pi*h^2, type='l', col='red')
legend("topleft", lwd=c(1.5,1), 
       col=c("black", "red"), 
       legend=c(expression(hat(K)), expression(pi*h^2)))
par(mar = c(5.1,4.1,4.1,2.1))
dev.off()

# Repeating k function for N2
kpts.2=khat(N2,Neuse.poly,h)

# k function plot
png("projects/results/sa_p3.png", width = 800, height = 600)
par(mar = c(5.1,4.1,1,1))
plot(h, kpts.2, type="l", lwd=1.5, ylab="K", xlab="Spatial Lag",
     main = "K-function vs. CSR curve | Weekday Group")
points(h, pi*h^2, type='l', col='red')
legend("topleft", lwd=c(1.5,1), 
       col=c("black", "red"), 
       legend=c(expression(hat(K)), expression(pi*h^2)))
par(mar = c(5.1,4.1,4.1,2.1))
dev.off()



# ______________________________________________________________________________
# Contour Intensity plots for Weekend and weekday data
# ______________________________________________________________________________
b=.0082 # bandwitdh of 0.5km
lam.nx=150
lam.ny=200

# intensity plot for N1
lam.est=kernel2d(N1, Neuse.poly, b, lam.nx, lam.ny)
png("projects/results/sa_p4.png", width = 800, height = 600)
image(lam.est$x,lam.est$y,lam.est$z,col=terrain.colors(100), asp=1,
      main="Map of Kernel Estimated Intensities on Weekends", xlab="Longitude", ylab="Latitude")
polymap(Neuse.poly,add=TRUE)
pointmap(N1,pch=20,cex=1,add=TRUE)
dev.off()

# intensity plot for N2
lam.est.2=kernel2d(N2, Neuse.poly, b, lam.nx, lam.ny)
png("projects/results/sa_p5.png", width = 800, height = 600)
image(lam.est.2$x,lam.est.2$y,lam.est.2$z,col=terrain.colors(100), asp=1,
      main="Map of Kernel Estimated Intensities on Weekdays", xlab="Longitude", ylab="Latitude")
polymap(Neuse.poly,add=TRUE)
pointmap(N2,pch=20,cex=1,add=TRUE)
SpectrumLegend("topleft", legend = c("Low","         ", "High"),
               palette = terrain.colors(100), lwd = 6,
               horiz = T, title="Intensity")
dev.off()

# ______________________________________________________________________________
# Comparisons of Point patterns between weekdays and weekends
# ______________________________________________________________________________
# estimating the log relative risk

# log of ratio of weekend intensity to weekday intensity
log.ratio <- log(lam.est$z / lam.est.2$z)

# estimation of the log ratio of normalized intensities
# assuming q0=q1=1
l2 <- log(n/n2)

# log relative risk estimate
log.risk <- log.ratio - l2

# plot of intensity comparison between weekend and weekday data
colpalette <- viridis_pal(option = "magma")
colpal_rev <- rev(colpalette(15))

png("projects/results/sa_p6.png", width = 800, height = 600)
image(lam.est$x,lam.est$y,log.risk,col=colpal_rev, asp=1,
      main="Map of Log Relative Risk", xlab="Longitude", ylab="Latitude",
      ylim = c(35.82, 35.95))
polymap(Neuse.poly, add=T)
SpectrumLegend("topleft", legend = c("Low","         ", "High"),
               palette = colpal_rev, lwd = 6,
               horiz = T, title="Relative Risk")
dev.off()

# ______________________________________________________________________________
# Monte Carlo test under a random labeling Hypothesis
# ______________________________________________________________________________
set.seed(608)

monte.sims <- 200

lrisk.sims <- list()

# point pattern for traffic crimes
tcrimes <- rbind(N1, N2)

for(i in 1:monte.sims) {
  print(i)
  
  samp <- sample(1:nrow(tcrimes), size = n, replace = F)
  
  # random samples from pool of traffic crimes
  weekend.sample <- tcrimes[samp,]
  weekday.sample <- tcrimes[-samp,]
  
  # intensity estimates
  intensity.weekend <- kernel2d(weekend.sample, Neuse.poly, b, lam.nx, lam.ny)
  intensity.weekday <- kernel2d(weekday.sample, Neuse.poly, b, lam.nx, lam.ny)
  
  lrisk.ratio <- log(intensity.weekend$z / intensity.weekday$z) - l2
  
  lrisk.sims[[i]] <- lrisk.ratio
}

lrisk.sims <- array(unlist(lrisk.sims), dim=c(lam.nx, lam.ny, monte.sims))

# tolerance intervals
alpha <- 0.05

lower.tol <- apply(lrisk.sims, c(1, 2), quantile,
                   probs = alpha/2, na.rm=T)
upper.tol <- apply(lrisk.sims, c(1, 2), quantile,
                   probs = 1-alpha/2, na.rm=T)

sig.vals <- matrix(0, lam.nx, lam.ny)
sig.vals[log.risk < lower.tol] <- -1
sig.vals[log.risk > upper.tol] <- 1

sig.width <- upper.tol-lower.tol

png("projects/results/sa_p7.png", width = 800, height = 600)
image(lam.est$x,lam.est$y,sig.vals,col=c("forestgreen", "white", "red3"), asp=1,
      main="Regions Exceeding 95% MC Intervals", xlab="Longitude", ylab="Latitude",
      ylim = c(35.82, 35.95))
polymap(Neuse.poly, add=T)
legend("topleft", legend = c("Higher Weekend Intensity", "Higher Weekday Intensity"),
       pch = 16, col = c("red3", "forestgreen"))
dev.off()

# spatial map of tolerance interval widths
image(lam.est$x,lam.est$y,sig.width,col=terrain.colors(100), asp=1,
      main="Sig Width", xlab="Longitude", ylab="Latitude",
      ylim = c(35.82, 35.95))
polymap(Neuse.poly, add=T)
legend_image <- as.raster(matrix(terrain.colors(100), ncol = 1))
rasterImage(legend_image, xleft = max(lam.est$x) + 0.1, ybottom = min(lam.est$y),
            xright = max(lam.est$x) + 0.3, ytop = max(lam.est$y))

