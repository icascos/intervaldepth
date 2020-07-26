source("https://raw.githubusercontent.com/icascos/intervaldepth/master/interval_depth.R")

# several interval depths w.r.t. a sample of 100 standard normal observations

par(mar=c(2.5,0,1,0))
plot(NA,xlim=c(-6,6),ylim=c(0,1),axes=F,ann=F)
axis(1, -5:5)

set.seed(123)
gaussian1 <- rnorm(100,0,1)
gaus1 <- cbind(gaussian1, rep(0.18, 100))
segments(-3.8, 0.18,  3.8, 0.18, col=176)
points(gaus1, pch=16, cex=0.56)
legend(3, 0.28, "N(0,1)", bty="n", cex=0.75)
zon<-m1m2(gaussian1)
segments(zon[1],0.16,zon[1],0.2,col=2)
segments(zon[2],0.16,zon[2],0.2,col=2)
legend(-6.5,0.28,paste("depth=",round(sm1m2depth(gaussian1,gaussian1),2),sep=""),bty="n",cex=0.8)

mu <- c(0,0,0,1)
sigma <- c(1,0.5,2,1)
height <- c(.36,.54,.72,.9)

for(i in 1:4) {
  gaussian2 <- rnorm(10,mu[i],sigma[i])
  gaus2 <- cbind(gaussian2, rep(height[i], 10))
  segments(-3.8, height[i], 3.8, height[i], col=176)
  gaus12 <- cbind(gaussian1, rep(height[i], 100))
  points(gaus12, pch=16, cex=0.56, col="honeydew4")
  points(gaus2, pch=17, cex=0.95, col=1)
  zon<-m1m2(gaussian2)
  segments(zon[1],height[i]-0.02,zon[1],height[i]+0.02,col=2)
  segments(zon[2],height[i]-0.02,zon[2],height[i]+0.02,col=2)
  legend(3, height[i]+.1,paste("N(",mu[i],",",sigma[i],")",sep=""), bty="n", cex=0.75)
  legend(-6.5,height[i]+.1,paste("depth=",round(sm1m2depth(gaussian2,gaussian1),2),sep=""),bty="n",cex=0.8)
}

title(main="Interval depths")

