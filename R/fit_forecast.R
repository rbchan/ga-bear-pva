## Female-only analysis

library(rjags)
library(lattice)
library(raster)

list.files()


## Load SCR data
load("../data/ga_bear_data_females.gzip")

## Inspect the loaded data objects
str(ch4d.f)
str(traps)

## Check that traps are in the same order in both objects
## Must be TRUE
all(rownames(traps)==dimnames(ch4d.f)[[3]])
all(rownames(oper)==rownames(traps))


## Format data

y <- ch4d.f ## 4D encounter histories
str(y)

y.dim <- dim(y)
nBears <- y.dim[1]
nTraps <- y.dim[2]
nWeeks <- y.dim[3]
nYears <- y.dim[4]
M <- 500 # Data aug. Maximum number of bears that could have been alive

all(rownames(prevcap.f)==rownames(y)) ## Must be TRUE

## Non-spatial, robust design format
y.nonsp.robust <- apply(y>0, c(1,3,4), any)*1L
str(y.nonsp.robust)

## Non-spatial, non-robust format
y.nonsp.nonrobust <- apply(y.nonsp.robust>0, c(1,3), any)*1L
str(y.nonsp.nonrobust)

## Spatial, non-robust format
y.sp.nonrobust <- apply(y>0, c(1,2,4), any)*1L
str(y.sp.nonrobust)


## First and last detections 
first.det <- apply(y.nonsp.nonrobust>0, 1, function(x) min(which(x)))
last.det <- apply(y.nonsp.nonrobust>0, 1, function(x) max(which(x)))
known.alive <- last.det-first.det+1
table(known.alive)


## Data augmentation
yaug <- array(0, c(M, nTraps, nWeeks, nYears))
yaug[1:nBears,,,] <- y

str(yaug)


## Known values of the alive/dead state matrix
zdata.aug <- matrix(NA, M, nYears)
for(i in 1:nBears) {
    zdata.aug[i,first.det[i]:last.det[i]] <- 1
}

## Visualize the observed information about lifetime
image(t(zdata.aug[1:nBears,]))

## Bears known to be alive in each year
colSums(zdata.aug[1:nBears,], na.rm=TRUE)



## Operational status of traps
nTrapsYr <- colSums(oper)
nTrapsYr
oper.traps <- matrix(0L, nTraps, nYears)
for(i in 1:nYears) {
    oper.traps[1:nTrapsYr[i],i] <- which(oper[,i]==1)
}



## State-space
x.range <- range(traps[,1])
y.range <- range(traps[,2])
buffer <- 4000
xlim0 <- x.range+c(-buffer, buffer)
ylim0 <- y.range+c(-buffer, buffer)

S360 <- raster("../data/state-space360.tif")

traps.in <- extract(S360, traps)

plot(S360)

## Convert raster to a data.frame
S360dat <- as.data.frame(S360, xy=TRUE)
names(S360dat) <- c("x", "y", "S360")
S360dat$x <- floor(S360dat$x)
S360dat$y <- floor(S360dat$y)
table(diff(S360dat$x))
table(diff(S360dat$y))

nPix360 <- nrow(S360dat)
dist360 <- matrix(NA, nPix360, nTraps)
for(j in 1:nTraps) {
    dist360[,j] <- sqrt((S360dat$x-traps[j,1])^2 +
                        (S360dat$y-traps[j,2])^2)
}

S360dat$pixIn <- (rowSums(dist360<buffer)>0)*1L
table(S360dat$pixIn)

dim.S360 <- dim(S360)

S360dat$habitat <- ifelse((S360dat$S360>0) &
                          (S360dat$pixIn>0),
                          1L, 0L)

## Visualize traps and state space
levelplot(habitat ~ x+y, S360dat, aspect='iso',
          panel=function(...) {
              panel.levelplot(...)
              lpoints(traps, pch=3)
          })

levelplot(S360+pixIn+habitat ~ x+y, S360dat, aspect='iso',
          panel=function(...) {
              panel.levelplot(...)
              lpoints(traps, pch=3)
          })

delta <- res(S360)[1]
xseq <- sort(unique(S360dat$x))
yseq <- sort(unique(S360dat$y), decreasing=TRUE)
xlim <- range(xseq)+c(-delta/2,delta/2)
ylim <- range(yseq)+c(-delta/2,delta/2)

dim.S360
range(ceiling((runif(1e5, xlim[1], xlim[2])-xlim[1])/delta))
range(ceiling((ylim[2]-runif(1e5, ylim[1], ylim[2]))/delta))


S <- matrix(S360dat$habitat, dim.S360[1], dim.S360[2],
            dimnames=list(yseq, xseq), byrow=TRUE)
str(S)

plot(raster(S))

sum(S)*delta^2 / 1e6 ## Area in km-sq
sum(S360dat$pixIn)*delta^2 / 1e6


## JAGS
## Models can take several hours to run

library(rjags)



### Put the data in a list
jd <- list(y=y, y.zero=matrix(0, M, nYears),
           z=zdata.aug,
           nBears=nBears, M=M,
           nTrapsYr=as.integer(nTrapsYr),
           nWeeks=nWeeks, nYears=nYears, nFuture=0,
           oper.traps=oper.traps,
           prevcap=prevcap.f+1,
           S=S, delta=delta,
           ones=matrix(1, M, nYears),
           x=traps, xlim=xlim, ylim=ylim)
jd$z=cbind(zdata.aug, matrix(NA, M, jd$nFuture))

str(jd)


#### Function to generate initial values
ji <- function() {
    nFuture <- jd$nFuture
    nBears <- jd$nBears
    M <- jd$M
    nYears <- jd$nYears
    nWeeks <- jd$nWeeks
    nFuture <- jd$nFuture
    traps <- jd$x
    nTraps <- nrow(traps)
    zi <- matrix(0, M, nYears+nFuture)
    ## Avoid N[t]=0, which causes problems with inital values of lambda
    if(nFuture>0)
        zi[(M-5):M,(nYears+1):(nYears+nFuture)] <- 1
    si1 <- traps[sample.int(nTraps, size=M, replace=TRUE,
                            prob=traps.in),]
    si.rec <- si <- array(NA, c(M, 2, nYears))
    si[,,1] <- si1
    si.rec[,,-1] <- si1
    for(i in 1:nBears) {
        zi[i,first.det[i]:last.det[i]] <- NA
        traplocs.i <- (rowSums(jd$y[i,,,])>0) * traps.in
        if(all(!traplocs.i))
            next
        si[i,,1] <- colMeans(traps[traplocs.i,,drop=FALSE])
        si.rec[i,,-1] <- colMeans(traps[traplocs.i,,drop=FALSE])
    }
    p0i <- array(NA, c(nWeeks, nYears, 2))
    p0i[,,2] <- matrix(runif(nWeeks*nYears), nWeeks, nYears)
    list(z=zi, s=si, s.rec=si.rec, p0=p0i,
         sigma=runif(1, 3000, 4000), behave=runif(1),
         mu.lphi=0, sig.lphi=.01,
         psi=c(runif(1), rep(NA, nYears+nFuture-1)),
         igamma0=runif(1, 0.2, 0.3), gamma1=runif(1, 0, 0.001))
}

str(ji())

## Visualize initial values of activity centers
plot(raster(S, xmn=xlim[1], ymn=ylim[1], xmx=xlim[2], ymx=ylim[2]))
points(ji()$s[,,1])
points(ji()$s.rec[,,2], pch=16, col=3)
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="blue")
#symbols(Sx, Sy, circles=rep(buffer, length(Sx)),
#        fg="blue", add=TRUE, inches=FALSE)


#### Parameters to monitor
jp <- c("p0", "sigma", "behave", "mu.lphi", "sig.lphi", 
        "igamma0", "gamma0", "gamma1", "phi", "gamma", "EN", 
        "N", "Ntotal", "survivors", "recruits", "lambda")


## Models in parallel with WAIC computation

library(parallel)

(nCores <- max(1, min(detectCores()-1, 4)))

cl1 <- makeCluster(4)

clusterExport(cl1, c("jd", "ji", "jp"

jm.out <- clusterEvalQ(cl1, {
    library(rjags)
    load.module("dic")
    load.module("glm")
    jm <- jags.model(file="dT-s0-rPDD2-phi0-pYB-sig0.jag",
                     data=jd, inits=ji,
                     n.chains=1, n.adapt=500)
    return(jm)
})


save(jm.out, file="jm21.gzip")



jc1.out <- clusterEvalQ(cl1, {
    jc1 <- coda.samples(model=jm,
                        variable.names=c(jp, "deviance"),
                        n.iter=1500)
    return(as.mcmc(jc1))
})


save(jc1.out, file="jc1-out21.gzip")


jc2.out <- clusterEvalQ(cl1, {
    jc2 <- coda.samples(model=jm,
                        variable.names=c(jp, "deviance", "z", "a","s"),
                        n.iter=5000, n.thin=1)
    return(as.mcmc(jc2))
})



save(jc2.out, file="jc2-out21.gzip")

rm(jc1.out, jc2.out)
gc()


## WAIC
js3.out <- clusterEvalQ(cl1, {
    js3 <- jags.samples(model=jm,
                        variable.names=c("deviance","WAIC"),
                        n.iter=1000, n.thin=1, type="mean")
    return(js3)
})


save(js3.out, file="js3-out21.gzip")


jc4.out <- clusterEvalQ(cl1, {
    jc4 <- coda.samples(model=jm,
                        variable.names=c(jp, "deviance", "z", "a","s"),
                        n.iter=5000, n.thin=1)
    return(as.mcmc(jc4))
})



save(jc4.out, file="jc4-out21.gzip")




## WAIC


js5.out <- clusterEvalQ(cl1, {
    js5 <- jags.samples(model=jm,
                        variable.names=c("deviance","WAIC"),
                        n.iter=1000, n.thin=1, type="mean")
    return(js5)
})


save(js5.out, file="js5-out21.gzip")




jm.out <- clusterEvalQ(cl1, {
    return(jm)
})


save(jm.out, file="jm21.gzip")



state.out <- clusterEvalQ(cl1, {
    return(jm$state())
})


save(state.out, file="jm-state21.gzip")























## Forecast

library(coda)

source("../forecast.R")

args(forecast)





## Model 21
## This is the model with density-dependent per-captita
## recruitment


library(coda)

source("../forecast.R")

load("mod21-no-sr/jc4-out21.gzip")
jc21f.4 <- as.mcmc.list(jc4.out)
jc21f.4 <- window(jc21f.4, thin=10)

niter(jc21f.4)*nchain(jc21f.4)

## debugonce(forecast)


test21 <- forecast(jc21f.4, nFuture=50, M=5000,
                   harvest=05, stochastic=TRUE,
                   maxRecruit=1,
                   NbarSims=1,
                   report=100)

Ntest <- apply(test21$z, c(2,3), sum)
matplot(Ntest, type="l")

EXtest <- apply(test21$Nbar==0, c(1,3), mean)

matplot(EXtest, type="l")

apply(EXtest, 1, mean)
apply(EXtest, 1, quantile, prob=0.975)




## No harvest
system.time({
    forc21.h0 <- forecast(jc21f.4, nFuture=50, M=15000,
                          harvest=0, stochastic=TRUE,
                          maxRecruit=1,
                          NbarSims=100,
                          report=100)
})

save(forc21.h0, file="forc21-h0.gzip")



## 5 additional females per year
system.time({
    forc21.h5 <- forecast(jc21f.4, nFuture=50, M=15000,
                          harvest=5, stochastic=TRUE,
                          maxRecruit=1,
                          NbarSims=100,
                          report=100)
})

save(forc21.h5, file="forc21-h5.gzip")



## 10 additional females per year
system.time({
    forc21.h10 <- forecast(jc21f.4, nFuture=50, M=15000,
                           harvest=10, stochastic=TRUE,
                           maxRecruit=1,
                           NbarSims=100,
                           report=100)
})

save(forc21.h10, file="forc21-h10.gzip")





## 15 additional females per year
system.time({
    forc21.h15 <- forecast(jc21f.4, nFuture=50, M=15000,
                           harvest=15, stochastic=TRUE,
                           maxRecruit=1,
                           NbarSims=100,
                           report=100)
})

save(forc21.h15, file="forc21-h15.gzip")




## 20 additional females per year
system.time({
    forc21.h20 <- forecast(jc21f.4, nFuture=50, M=15000,
                           harvest=20, stochastic=TRUE,
                           maxRecruit=1,
                           NbarSims=100,
                           report=100)
})

save(forc21.h20, file="forc21-h20.gzip")





## 25 additional females per year
system.time({
    forc21.h25 <- forecast(jc21f.4, nFuture=50, M=15000,
                           harvest=25, stochastic=TRUE,
                           maxRecruit=1,
                           NbarSims=100,
                           report=100)
})

save(forc21.h25, file="forc21-h25.gzip")




load("forc21-h0.gzip")
load("forc21-h5.gzip")
load("forc21-h10.gzip")
load("forc21-h15.gzip")
load("forc21-h20.gzip")
load("forc21-h25.gzip")



## Extinction risk curves for each posterior draw
EX21.h0 <- apply(forc21.h0$Nbar==0, c(1, 3), mean)
EX21.h5 <- apply(forc21.h5$Nbar==0, c(1, 3), mean)
EX21.h10 <- apply(forc21.h10$Nbar==0, c(1, 3), mean)
EX21.h15 <- apply(forc21.h15$Nbar==0, c(1, 3), mean)
EX21.h20 <- apply(forc21.h20$Nbar==0, c(1, 3), mean)
EX21.h25 <- apply(forc21.h25$Nbar==0, c(1, 3), mean)

matplot(EX21.h25[,1:100])

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
EN21.h0 <- apply(forc21.h0$Nbar, c(1, 3), median)
EN21.h5 <- apply(forc21.h5$Nbar, c(1, 3), median)
EN21.h10 <- apply(forc21.h10$Nbar, c(1, 3), median)
EN21.h15 <- apply(forc21.h15$Nbar, c(1, 3), median)
EN21.h20 <- apply(forc21.h20$Nbar, c(1, 3), median)
EN21.h25 <- apply(forc21.h25$Nbar, c(1, 3), median)

matplot(EN21.h0, type="l")
matplot(EN21.h15, type="l")
matplot(EN21.h25, type="l")



quantile(EN21.h0[51,], c(0.025, 0.5, 0.975))
quantile(EN21.h5[51,], c(0.025, 0.5, 0.975))
quantile(EN21.h10[51,], c(0.025, 0.5, 0.975))
quantile(EN21.h15[51,], c(0.025, 0.5, 0.975))
quantile(EN21.h20[51,], c(0.025, 0.5, 0.975))



fyrs <- 2017:(2017+50)




## Harvest rates (crude way of computing this)
5/rowMeans(EN21.h5)
10/rowMeans(EN21.h10)
15/rowMeans(EN21.h15)
20/rowMeans(EN21.h20)





png("fig/Erisk21.png", width=7, height=5, units="in", res=400)
par(mai=c(0.85, 0.85, 0.1, 0.1))
plot(fyrs, EX21.h0.mean, type="o", col="red", pch=16, cex=0.7, lwd=1.5, ylim=c(0,1),
     xlab="Year", ylab="Probability of metapopulation extinction", cex.lab=1.1)
polygon(c(fyrs, rev(fyrs)), c(EX21.h25.low, rev(EX21.h25.upp)),
        border=NA, col=rgb(1,0.15,0,0.6))
polygon(c(fyrs, rev(fyrs)), c(EX21.h20.low, rev(EX21.h20.upp)),
        border=NA, col=rgb(1,0.3,0,0.6))
polygon(c(fyrs, rev(fyrs)), c(EX21.h15.low, rev(EX21.h15.upp)),
        border=NA, col=rgb(1,0.45,0,0.6))
polygon(c(fyrs, rev(fyrs)), c(EX21.h10.low, rev(EX21.h10.upp)),
        border=NA, col=rgb(1,0.6,0,1))
polygon(c(fyrs, rev(fyrs)), c(EX21.h5.low, rev(EX21.h5.upp)),
        border=NA, col=rgb(1,0.75,0,1))
polygon(c(fyrs, rev(fyrs)), c(EX21.h0.low, rev(EX21.h0.upp)),
        border=NA, col=rgb(1,0.9,0,1))
lines(fyrs, EX21.h25.mean, type="o", col=rgb(1,0.15,0), #"red3",
      pch=15, cex=0.7, lty=1, lwd=1.5)
lines(fyrs, EX21.h20.mean, type="o", col=rgb(1,0.3,0), #"red3",
      pch=15, cex=0.7, lty=1, lwd=1.5)
lines(fyrs, EX21.h15.mean, type="o", col=rgb(1,0.45,0), #"red3",
      pch=15, cex=0.7, lty=1, lwd=1.5)
lines(fyrs, EX21.h10.mean, type="o", col=rgb(1,0.6,0), #"red2",
      pch=18, cex=0.7, lty=1, lwd=1.5)
lines(fyrs, EX21.h5.mean, type="o", col=rgb(1,0.75,0), #"red1",
      pch=17, cex=0.7, lty=1, lwd=1.5)
lines(fyrs, EX21.h0.mean, type="o", col=rgb(1,0.9,0), #"red1",
      pch=16, cex=0.7, lty=1, lwd=1.5)
legend(fyrs[1], 0.4, #paste(c(0, 1, 3, 6), "Status quo"),
       rev(c("5 additional females harvested", "10 additional females harvested",
             "15 additional females harvested", "20 additional females harvested",
             "25 additional females harvested")),
##       col=c("red", "red1", "red2", "red3"),
       col=c(rgb(1,0.15,0), rgb(1,0.3,0), rgb(1,0.45,0), rgb(1,0.6,0),
             rgb(1,0.75,0), rgb(1,0.9,0)),
       lty=1,
       pt.cex=0.9, pch=rev(c(16,17,18,15,19,20)))
dev.off()
system("open fig/Erisk21.png")





png("fig/Erisk21-noCI.png", width=7, height=6, units="in", res=400)
par(mai=c(0.85, 0.85, 0.1, 0.1))
plot(fyrs, EX21.h0.mean, type="o", col="red", pch=16, cex=0.7, lwd=1.5, ylim=c(0,1),
     xlab="Year", ylab="Extinction probability", cex.lab=1.1)
lines(fyrs, EX21.h25.mean, type="o", col=rgb(1,0.15,0), #"red3",
      pch=13, cex=0.7, lty=1, lwd=1.5)
lines(fyrs, EX21.h20.mean, type="o", col=rgb(1,0.3,0), #"red3",
      pch=14, cex=0.7, lty=1, lwd=1.5)
lines(fyrs, EX21.h15.mean, type="o", col=rgb(1,0.45,0), #"red3",
      pch=15, cex=0.7, lty=1, lwd=1.5)
lines(fyrs, EX21.h10.mean, type="o", col=rgb(1,0.6,0), #"red2",
      pch=18, cex=0.7, lty=1, lwd=1.5)
lines(fyrs, EX21.h5.mean, type="o", col=rgb(1,0.75,0), #"red1",
      pch=17, cex=0.7, lty=1, lwd=1.5)
lines(fyrs, EX21.h0.mean, type="o", col=rgb(1,0.9,0), #"red1",
      pch=16, cex=0.7, lty=1, lwd=1.5)
legend(fyrs[1], 1, #paste(c(0, 1, 3, 6), "Status quo"),
       rev(c("Status quo",
             "5 additional females harvested", "10 additional females harvested",
             "15 additional females harvested", "20 additional females harvested",
             "25 additional females harvested")),
##       col=c("red", "red1", "red2", "red3"),
       col=c(rgb(1,0.15,0), rgb(1,0.3,0), rgb(1,0.45,0), rgb(1,0.6,0),
             rgb(1,0.75,0), rgb(1,0.9,0)),
       lty=1,
       pt.cex=0.9, pch=rev(c(16,17,18,15,14,13)))
dev.off()
system("open fig/Erisk21-noCI.png")





## Stochastic post-repro harvest of 30 females/yr
png("fig/forc21-h0.png", width=8, height=6, units="in", res=400)
##pdf("fig/proj20-h30s.pdf", width=8, height=6)
##par(mai=c(0.9, 0.9, 0.9, 0.9))
matplot(fyrs, EN21.h0, type="l", xlab="Year",
        ylab="Female abundance", ylim=c(0, 1000),
        cex.lab=1.5, col=gray(0.8),
        lwd=0.5)
polygon(c(fyrs, rev(fyrs)),
        c(apply(EN21.h0, 1, quantile, prob=0.025),
          rev(apply(EN21.h0, 1, quantile, prob=0.975))),
        border=NA, col=rgb(1,1,0,0.6))
lines(fyrs, rowMeans(EN21.h0), lwd=4, col="yellow")
dev.off()
system("open fig/forc21-h0.png")





fyrs <- 2017:(2017+50)
iyrs <- 1:51 ## Reduce the length of this sequence to shorten the
             ## timehorizon in the figure


pclr <- rgb(0,0,1,0.6)
lclr <- rgb(0,0,1)

png("fig/forc21-h0-h20.png", width=6, height=9, units="in", res=400)
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
plot(fyrs[iyrs], EX21.h0.mean[iyrs], type="l", col="red", lwd=3, ylim=c(0,0.5),
     xlab="", ylab="Extinction risk", cex.lab=1.3)
## polygon(c(fyrs[iyrs], rev(fyrs[iyrs])), c(EX21.h0.low[iyrs], rev(EX21.h0.upp[iyrs])),
##         border=NA, col=rgb(1,0.15,0,0.6))
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
plot(fyrs[iyrs], EX21.h5.mean[iyrs], type="l", col="red", lwd=3, ylim=c(0,0.5),
     xlab="", ylab="Extinction risk", cex.lab=1.3)
## polygon(c(fyrs[iyrs], rev(fyrs[iyrs])), c(EX21.h5.low[iyrs], rev(EX21.h5.upp[iyrs])),
##         border=NA, col=rgb(1,0.15,0,0.6))
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
plot(fyrs[iyrs], EX21.h10.mean[iyrs], type="l", col="red", lwd=3, ylim=c(0,0.5),
     xlab="", ylab="Extinction risk", cex.lab=1.3)
## polygon(c(fyrs[iyrs], rev(fyrs[iyrs])), c(EX21.h10.low[iyrs], rev(EX21.h10.upp[iyrs])),
##         border=NA, col=rgb(1,0.15,0,0.6))
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
plot(fyrs[iyrs], EX21.h15.mean[iyrs], type="l", col="red", lwd=3, ylim=c(0,0.5),
     xlab="", ylab="Extinction risk", cex.lab=1.3)
## polygon(c(fyrs[iyrs], rev(fyrs[iyrs])), c(EX21.h15.low[iyrs], rev(EX21.h15.upp[iyrs])),
##         border=NA, col=rgb(1,0.15,0,0.6))
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
plot(fyrs[iyrs], EX21.h20.mean[iyrs], type="l", col="red", lwd=3, ylim=c(0,0.5),
     xlab="Year", ylab="Extinction risk", cex.lab=1.3)
## polygon(c(fyrs[iyrs], rev(fyrs[iyrs])), c(EX21.h20.low[iyrs], rev(EX21.h20.upp[iyrs])),
##         border=NA, col=rgb(1,0.15,0,0.6))
mtext("20 additional females harvested each year", side=3, line=0.4,
      outer=F, cex=1, at=c(2005,1))
text(2047, 0.45, "Extinction risk at year 50 = 37.0%")
dev.off()
system("open fig/forc21-h0-h20.png")



ls()

load("mod21-no-sr/jc2-out21.gzip")
jc21.2 <- as.mcmc.list(jc2.out)

vn21 <- varnames(jc21.2)

vn21sub <- c(vn21[-c(grep("p0", vn21), grep("z\\[", vn21),
                     grep("s\\[", vn21), grep("a\\[", vn21))],
             "p0[1,1,1]", "p0[1,1,2]", "p0[1,2,1]", "p0[1,2,2]",
             "p0[1,3,1]", "p0[1,3,2]", "p0[1,4,1]", "p0[1,4,2]",
             "p0[1,5,1]", "p0[1,5,2]",
             paste("gamma[", 1:4, "]", sep=""),
             paste("survivors[", 1:4, "]", sep=""),
             paste("recruits[", 1:4, "]", sep=""),
             paste("lambda[", 1:4, "]", sep=""), "sigma")

##vn21sub <- vn21sub[-grep("deviance", vn21sub)]

mc21.2sub <- as.data.frame(as.matrix(jc21.2[,vn21sub]))
summary(mc21.2sub)

## Area of the state-space
## sum(jdata15$S)*jdata15$delta^2/1e6
SSarea <- 622.5984 ## sq-km



## Expected density (bears per 100 sq-km)
mc21.2sub$"ED[1]" <- mc21.2sub$"EN[1]"/SSarea*100
mc21.2sub$"ED[2]" <- mc21.2sub$"EN[2]"/SSarea*100
mc21.2sub$"ED[3]" <- mc21.2sub$"EN[3]"/SSarea*100
mc21.2sub$"ED[4]" <- mc21.2sub$"EN[4]"/SSarea*100
mc21.2sub$"ED[5]" <- mc21.2sub$"EN[5]"/SSarea*100

## ## Baseline cap probs
## mc21.2sub$p0pre <- mc21.2sub$"p0[1,1,1]"
## mc21.2sub$p0post <- mc21.2sub$"p0[1,1,1]" / mc21.2sub$behave

summary(mc21.2sub)


mc21.2sub <- data.matrix(mc21.2sub)



sstats21.2.0 <- cbind(Mean=colMeans(mc21.2sub),
                      SD=apply(mc21.2sub, 2, sd),
                      LowerCI=apply(mc21.2sub, 2, quantile, prob=0.025),
                      Median=apply(mc21.2sub, 2, quantile, prob=0.5),
                      UpperCI=apply(mc21.2sub, 2, quantile, prob=0.975))

vn21out <- c(paste("N[", 1:5, "]", sep=""),
             paste("ED[", 1:5, "]", sep=""),
             "Ntotal",
             "phi",
             paste("survivors[", 1:4, "]", sep=""),
             "igamma0", "gamma1",
             paste("gamma[", 1:4, "]", sep=""),
             paste("recruits[", 1:4, "]", sep=""),
             paste("lambda[", 1:4, "]", sep=""),
             ##             "p0pre", "p0post",
             grep("p0", rownames(sstats21.2.0), value=TRUE),
             "sigma")




sstats21.2 <- sstats21.2.0[vn21out,]


sstats21.2

write.table(format(sstats21.2, digits=1, sci=FALSE),
            quote=FALSE, sep="\t")

write.table(format(sstats21.2, digits=1, sci=FALSE),
            file="sstats21-2.txt",
            quote=FALSE, sep="\t")


write.table(round(sstats21.2, digits=2),
            quote=FALSE, sep="\t")

write.table(round(sstats21.2, digits=2),
            file="sstats21-2.txt",
            quote=FALSE, sep="\t")






png("fig/N-21.png", width=8, height=6, units="in", res=400)
par(mai=c(0.9, 0.9, 0.9, 0.9))
yrs <- 2012:2016
plot(yrs, sstats21.2.0[paste("N[", 1:5, "]", sep=""),"Median"],
     type="b", ylim=c(70, 210), pch=16,
     xlab="Year", ylab="Female abundance", cex.lab=1.5)
axis(4, at=seq(80, 200, 20),
     label=round(seq(80, 200, 20)/SSarea*100, 1))
mtext(expression(paste("Density (females / 100 ", km^2, ")")),
      side=4, line=-1.5,
      outer=TRUE, cex=1.5)
arrows(yrs,
       sstats21.2.0[paste("N[", 1:5, "]", sep=""),"LowerCI"],
       yrs,
       sstats21.2.0[paste("N[", 1:5, "]", sep=""),"UpperCI"],
       angle=90, code=3, length=0.05)
dev.off()
system("open fig/N-21.png")



## Density dependence

gamma01 <- as.matrix(jc21.2[,c("igamma0","gamma1")])[seq(1, 20000, by=10),]

str(gamma01)

Nx <- seq(0, 200)

DDpost <- apply(gamma01, 1, function(x) x[1]*exp(-x[2]*Nx^2))

DDpost.mean <- apply(DDpost, 1, mean)
DDpost.low <- apply(DDpost, 1, quantile, prob=0.025)
DDpost.upp <- apply(DDpost, 1, quantile, prob=0.975)


##png("fig/DDR-21.png", width=6, height=5, units="in", res=500)
pdf("fig/DDR-21.pdf", width=6, height=5)
par(mai=c(0.8, 0.8, 0.1, 0.1))
plot(sstats21.2.0[paste("N[", 1:4, "]", sep=""),"Median"],
     sstats21.2.0[paste("gamma[", 1:4, "]", sep=""),"Median"],
     xlim=c(0, 200), xaxt="n",
     ylim=c(0, 1),
     xlab=expression(paste("Density (females / 100 ", km^2, ")")),
##     xlab="Female abundance",
     ylab="Per-capita recruitment")
##matlines(DDpost, type="l", col=gray(0.8))
points(sstats21.2.0[paste("N[", 1:4, "]", sep=""),"Median"],
       sstats21.2.0[paste("gamma[", 1:4, "]", sep=""),"Median"], pch=16)
arrows(sstats21.2.0[paste("N[", 1:4, "]", sep=""),"Median"],
       sstats21.2.0[paste("gamma[", 1:4, "]", sep=""),"LowerCI"],
       sstats21.2.0[paste("N[", 1:4, "]", sep=""),"Median"],
       sstats21.2.0[paste("gamma[", 1:4, "]", sep=""),"UpperCI"],
       angle=90, code=3, length=0.05)
lines(Nx, DDpost.mean)
lines(Nx, DDpost.low, lty=2)
lines(Nx, DDpost.upp, lty=2)
Dx <- seq(0, 30, 10)
##axis(1, at=seq(0, 200, 50), labels=round(seq(0, 200, 50)/SSarea*100, 1))
axis(1, at=Dx*SSarea/100, labels=Dx)
dev.off()
##system("open fig/DDR-21.png")
system("open fig/DDR-21.pdf")




