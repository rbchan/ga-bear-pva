## Forecast Central GA black bear dynamics and viability

library(coda)

source("forecast_fn.R")

args(forecast)


## Load posterior samples
load("jc4.gzip")
jc4 <- as.mcmc.list(jc4.out)
ps <- window(jc4, thin=10)

niter(ps)*nchain(ps)




## debugonce(forecast)


## No harvest
system.time({
    forc.h0 <- forecast(ps, nFuture=50, M=15000,
                        harvest=0, stochastic=TRUE,
                        NbarSims=100,
                        report=100)
})

save(forc.h0, file="forc-h0.gzip")



## 5 additional females per year
system.time({
    forc.h5 <- forecast(ps, nFuture=50, M=15000,
                        harvest=5, stochastic=TRUE,
                        NbarSims=100,
                        report=100)
})

save(forc.h5, file="forc-h5.gzip")



## 10 additional females per year
system.time({
    forc.h10 <- forecast(ps, nFuture=50, M=15000,
                         harvest=10, stochastic=TRUE,
                         NbarSims=100,
                         report=100)
})

save(forc.h10, file="forc-h10.gzip")





## 15 additional females per year
system.time({
    forc.h15 <- forecast(ps, nFuture=50, M=15000,
                         harvest=15, stochastic=TRUE,
                         NbarSims=100,
                         report=100)
})

save(forc.h15, file="forc-h15.gzip")




## 20 additional females per year
system.time({
    forc.h20 <- forecast(ps, nFuture=50, M=15000,
                         harvest=20, stochastic=TRUE,
                         NbarSims=100,
                         report=100)
})

save(forc.h20, file="forc-h20.gzip")





## 25 additional females per year
system.time({
    forc.h25 <- forecast(ps, nFuture=50, M=15000,
                         harvest=25, stochastic=TRUE,
                         NbarSims=100,
                         report=100)
})

save(forc.h25, file="forc-h25.gzip")



## Load saved results (if the forecasting above has already been completed)
load("forc-h0.gzip")
load("forc-h5.gzip")
load("forc-h10.gzip")
load("forc-h15.gzip")
load("forc-h20.gzip")
load("forc-h25.gzip")



## Extinction risk curves for each posterior draw
EX21.h0 <- apply(forc.h0$Nbar==0, c(1, 3), mean)
EX21.h5 <- apply(forc.h5$Nbar==0, c(1, 3), mean)
EX21.h10 <- apply(forc.h10$Nbar==0, c(1, 3), mean)
EX21.h15 <- apply(forc.h15$Nbar==0, c(1, 3), mean)
EX21.h20 <- apply(forc.h20$Nbar==0, c(1, 3), mean)
EX21.h25 <- apply(forc.h25$Nbar==0, c(1, 3), mean)

matplot(EX21.h25[,1:100])

## Posterior means and 95% CIs
EX21.h0.mean <- rowMeans(EX21.h0)
EX21.h0.low <- apply(EX21.h0, 1, quantile, prob=0.025)
EX21.h0.upp <- apply(EX21.h0, 1, quantile, prob=0.975)

EX21.h5.mean <- rowMeans(EX21.h5)
EX21.h5.low <- apply(EX21.h5, 1, quantile, prob=0.025)
EX21.h5.upp <- apply(EX21.h5, 1, quantile, prob=0.975)

EX21.h10.mean <- rowMeans(EX21.h10)
EX21.h10.low <- apply(EX21.h10, 1, quantile, prob=0.025)
EX21.h10.upp <- apply(EX21.h10, 1, quantile, prob=0.975)

EX21.h15.mean <- rowMeans(EX21.h15)
EX21.h15.low <- apply(EX21.h15, 1, quantile, prob=0.025)
EX21.h15.upp <- apply(EX21.h15, 1, quantile, prob=0.975)

EX21.h20.mean <- rowMeans(EX21.h20)
EX21.h20.low <- apply(EX21.h20, 1, quantile, prob=0.025)
EX21.h20.upp <- apply(EX21.h20, 1, quantile, prob=0.975)

EX21.h25.mean <- rowMeans(EX21.h25)
EX21.h25.low <- apply(EX21.h25, 1, quantile, prob=0.025)
EX21.h25.upp <- apply(EX21.h25, 1, quantile, prob=0.975)


plot(EX.h25.mean, ylim=0:1)
lines(EX.h25.upp)




## Annual abundance for each posterior draw, averaged
## over stochastic harvest
EN21.h0 <- apply(forc.h0$Nbar, c(1, 3), median)
EN21.h5 <- apply(forc.h5$Nbar, c(1, 3), median)
EN21.h10 <- apply(forc.h10$Nbar, c(1, 3), median)
EN21.h15 <- apply(forc.h15$Nbar, c(1, 3), median)
EN21.h20 <- apply(forc.h20$Nbar, c(1, 3), median)
EN21.h25 <- apply(forc.h25$Nbar, c(1, 3), median)

matplot(EN21.h0, type="l")
matplot(EN21.h15, type="l")
matplot(EN21.h25, type="l")



quantile(EN21.h0[51,], c(0.025, 0.5, 0.975))
quantile(EN21.h5[51,], c(0.025, 0.5, 0.975))
quantile(EN21.h10[51,], c(0.025, 0.5, 0.975))
quantile(EN21.h15[51,], c(0.025, 0.5, 0.975))
quantile(EN21.h20[51,], c(0.025, 0.5, 0.975))



## Harvest rates 
5/rowMeans(EN21.h5)
10/rowMeans(EN21.h10)
15/rowMeans(EN21.h15)
20/rowMeans(EN21.h20)




## Years
fyrs <- 2017:(2017+50)
iyrs <- 1:51 ## Reduce the length of this sequence to shorten the
             ## timehorizon in the figure

## Colors
pclr <- rgb(0,0,1,0.6)
lclr <- rgb(0,0,1)

## Figure 6 in the paper
png("../figs/fig-6_forecasts_h0-h20.png", width=6, height=9,
    units="in", res=400)
par(mfrow=c(5,2), mai=c(0.6, 0.6, 0.2, 0.1))
## h0
matplot(fyrs[iyrs], EN21.h0[iyrs,], type="l", xlab="",
        ylab="Female abundance", ylim=c(0, 600),
        col=gray(0.8), cex.lab=1.3,
        lwd=0.5)
polygon(c(fyrs[iyrs], rev(fyrs[iyrs])),
        c(apply(EN21.h0[iyrs,], 1, quantile, prob=0.025),
          rev(apply(EN21.h0[iyrs,], 1, quantile, prob=0.975))),
        border=NA, col=pclr) #col=rgb(1,1,0,0.6))
lines(fyrs[iyrs], rowMeans(EN21.h0)[iyrs], lwd=3, col=lclr)
plot(fyrs[iyrs], EX21.h0.mean[iyrs], type="l", col="red", lwd=3,
     ylim=c(0,0.5),
     xlab="", ylab="Extinction risk", cex.lab=1.3)
mtext("Status quo", side=3, line=0.4, outer=F, cex=1, at=c(2005,1))
text(2047, 0.45, "Extinction risk at year 50 < 0.001%")
## h5
matplot(fyrs[iyrs], EN21.h5[iyrs,], type="l", xlab="",
        ylab="Female abundance", ylim=c(0, 600),
        col=gray(0.8), cex.lab=1.3,
        lwd=0.5)
polygon(c(fyrs[iyrs], rev(fyrs[iyrs])),
        c(apply(EN21.h5[iyrs,], 1, quantile, prob=0.025),
          rev(apply(EN21.h5[iyrs,], 1, quantile, prob=0.975))),
        border=NA, col=pclr) #col=rgb(1,1,0,0.6))
lines(fyrs[iyrs], rowMeans(EN21.h5)[iyrs], lwd=3, col=lclr)
plot(fyrs[iyrs], EX21.h5.mean[iyrs], type="l", col="red", lwd=3,
     ylim=c(0,0.5),
     xlab="", ylab="Extinction risk", cex.lab=1.3)
mtext("5 additional females harvested each year", side=3, line=0.4,
      outer=F, cex=1, at=c(2005,1))
text(2047, 0.45, "Extinction risk at year 50 = 0.13%")
## h10
matplot(fyrs[iyrs], EN21.h10[iyrs,], type="l", xlab="",
        ylab="Female abundance", ylim=c(0, 600),
        col=gray(0.8), cex.lab=1.3,
        lwd=0.5)
polygon(c(fyrs[iyrs], rev(fyrs[iyrs])),
        c(apply(EN21.h10[iyrs,], 1, quantile, prob=0.025),
          rev(apply(EN21.h10[iyrs,], 1, quantile, prob=0.975))),
        border=NA, col=pclr) #col=rgb(1,1,0,0.6))
lines(fyrs[iyrs], rowMeans(EN21.h10)[iyrs], lwd=3, col=lclr)
plot(fyrs[iyrs], EX21.h10.mean[iyrs], type="l", col="red", lwd=3,
     ylim=c(0,0.5),
     xlab="", ylab="Extinction risk", cex.lab=1.3)
mtext("10 additional females harvested each year", side=3, line=0.4,
      outer=F, cex=1, at=c(2005,1))
text(2047, 0.45, "Extinction risk at year 50 = 1.15%")
## h15
matplot(fyrs[iyrs], EN21.h15[iyrs,], type="l", xlab="",
        ylab="Female abundance", ylim=c(0, 600),
        col=gray(0.8), cex.lab=1.3,
        lwd=0.5)
polygon(c(fyrs[iyrs], rev(fyrs[iyrs])),
        c(apply(EN21.h15[iyrs,], 1, quantile, prob=0.025),
          rev(apply(EN21.h15[iyrs,], 1, quantile, prob=0.975))),
        border=NA, col=pclr) #col=rgb(1,1,0,0.6))
lines(fyrs[iyrs], rowMeans(EN21.h15)[iyrs], lwd=3, col=lclr)
plot(fyrs[iyrs], EX21.h15.mean[iyrs], type="l", col="red", lwd=3,
     ylim=c(0,0.5),
     xlab="", ylab="Extinction risk", cex.lab=1.3)
mtext("15 additional females harvested each year", side=3, line=0.4,
      outer=F, cex=1, at=c(2005,1))
text(2047, 0.45, "Extinction risk at year 50 = 9.0%")
## h20
matplot(fyrs[iyrs], EN21.h20[iyrs,], type="l", xlab="Year",
        ylab="Female abundance", ylim=c(0, 600),
        col=gray(0.8), cex.lab=1.3,
        lwd=0.5)
polygon(c(fyrs[iyrs], rev(fyrs[iyrs])),
        c(apply(EN21.h20[iyrs,], 1, quantile, prob=0.025),
          rev(apply(EN21.h20[iyrs,], 1, quantile, prob=0.975))),
        border=NA, col=pclr) #col=rgb(1,1,0,0.6))
lines(fyrs[iyrs], rowMeans(EN21.h20)[iyrs], lwd=3, col=lclr)
plot(fyrs[iyrs], EX21.h20.mean[iyrs], type="l", col="red", lwd=3,
     ylim=c(0,0.5),
     xlab="Year", ylab="Extinction risk", cex.lab=1.3)
mtext("20 additional females harvested each year", side=3, line=0.4,
      outer=F, cex=1, at=c(2005,1))
text(2047, 0.45, "Extinction risk at year 50 = 37.0%")
dev.off()
system("open ../figs/fig-6_forecasts_h0-h20.png")



