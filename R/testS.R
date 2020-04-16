## This script assesses the sensitivity of the results to the
## specification of the state-space.
## It takes at least a day or two to run everything in this script.


library(rjags)
library(lattice)
library(raster)

list.files()

### ---------------------- Prepare the data ------------------------

## Load SCR data
load("../data/ga_bear_data_females.gzip")

## Inspect loaded data objects
str(ch4d.f)
str(traps)

## Check that traps are in the same order in both objects
## Must be TRUE
all(rownames(traps)==dimnames(ch4d.f)[[3]])
all(rownames(oper)==rownames(traps))


## Format data

y <- ch4d.f
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

## Determine the first and last years when each bear was detected
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
buffer2500 <- 2500
xlim0 <- x.range+c(-buffer, buffer)
ylim0 <- y.range+c(-buffer, buffer)

S360 <- raster("../data/state-space360.tif")

## Determine which traps were in the state-space
## Used later on to generate initial values of the activity centers
traps.in <- extract(S360, traps)


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


S360dat$Original <- ifelse((S360dat$S360>0) &
                           (S360dat$pixIn>0),
                           1L, 0L)

hab360dat <- subset(S360dat, Original==1)



delta <- res(S360)[1]
xseq <- sort(unique(S360dat$x))
yseq <- sort(unique(S360dat$y), decreasing=TRUE)
xlim <- range(xseq)+c(-delta/2,delta/2)
ylim <- range(yseq)+c(-delta/2,delta/2)

dim.S360
range(ceiling((runif(1e5, xlim[1], xlim[2])-xlim[1])/delta))
range(ceiling((ylim[2]-runif(1e5, ylim[1], ylim[2]))/delta))


## Original state-space
S1 <- matrix(S360dat$Original, dim.S360[1], dim.S360[2],
             dimnames=list(yseq, xseq), byrow=TRUE)
str(S1)


## State-space with holes filled in
S2 <- S1
for(i in 2:(nrow(S2)-1)) {
    for(j in 2:(ncol(S2)-1)) {
        if(S1[i,j]>0)
            next
        any.north <- any(S1[1:(i-1),j]>0)
        any.south <- any(S1[(i+1):nrow(S2),j]>0)
        any.west <- any(S1[i,1:(j-1)]>0)
        any.east <- any(S1[i,(j+1):ncol(S2)]>0)
        if(any.north & any.south & any.west & any.east)
            S2[i,j] <- 1L
    }
}


S12 <- stack(raster(S1), raster(S2))
names(S12) <- c("Original", "Filled")

plot(S12, axes=FALSE)
            

S12df <- as.data.frame(S12, xy=TRUE)

S360dat$Filled <- as.numeric(t(S2))

## Two specifications of state-space
png("../figs/state-spaces.png", width=7, height=5,
     units="in", res=400)
levelplot(Original+Filled ~ x+y, S360dat, aspect='iso',
          panel=function(...) {
              panel.levelplot(...)
              lpoints(traps, pch=3, col="blue", cex=0.5)
          col.regions=c("gray90", "seagreen3"),
          colorkey=FALSE, 
          layout=c(2,1), xlab="UTM East (m)", ylab="UTM North (m)")
dev.off()
system("open ../figs/state-spaces.png")

ls()

### --------------- Run the model with each state-space -------------



## JAGS
## These models can take >2 hours to run

library(rjags)


### Put the data in a list
jd1 <- list(y=y, y.zero=matrix(0, M, nYears),
           z=zdata.aug,
           nBears=nBears, M=M,
           nTrapsYr=as.integer(nTrapsYr),
           nWeeks=nWeeks, nYears=nYears, nFuture=0,
           oper.traps=oper.traps,
           prevcap=prevcap.f+1,
           S=S1, delta=delta,
           ones=matrix(1, M, nYears),
           x=traps, xlim=xlim, ylim=ylim)
jd1$z=cbind(zdata.aug, matrix(NA, M, jd1$nFuture))

str(jd1)


#### Function to generate initial values
ji <- function() {
    nFuture <- jd1$nFuture
    nBears <- jd1$nBears
    M <- jd1$M
    nYears <- jd1$nYears
    nWeeks <- jd1$nWeeks
    nFuture <- jd1$nFuture
    traps <- jd1$x
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
        traplocs.i <- (rowSums(jd1$y[i,,,])>0) * traps.in
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
plot(raster(S1, xmn=xlim[1], ymn=ylim[1], xmx=xlim[2], ymx=ylim[2]))
points(ji()$s[,,1])
points(ji()$s.rec[,,2], pch=16, col=3)
rect(xlim[1], ylim[1], xlim[2], ylim[2], border="blue")#
## symbols(Sx, Sy, circles=rep(buffer, length(Sx)),
##         fg="blue", add=TRUE, inches=FALSE)


#### Parameters to monitor
jp <- c("p0", "sigma", "behave", "mu.lphi", "sig.lphi", 
        "igamma0", "gamma0", "gamma1", "phi", "gamma", "EN", 
        "N", "Ntotal", "survivors", "recruits", "lambda")


## Models in parallel 

library(parallel)

(nCores <- max(1, min(detectCores()-1, 4)))

cl1 <- makeCluster(nCores)

clusterExport(cl1, c("jd1", "ji", "jp",
                     "traps.in", "first.det", "last.det"))

## Compile and adapt
system.time({
    jm.out <- clusterEvalQ(cl1, {
        library(rjags)
        load.module("dic")
        load.module("glm")
        jm <- jags.model(file="dT-s0-rPDD2-phi0-pYB-sig0.jag",
                         data=jd1, inits=ji,
                         n.chains=1, n.adapt=500)
        return(jm)
    })
}) ## 5904s on 3.1 GHz processor




## Draw posterior samples
system.time({
    jc1.out <- clusterEvalQ(cl1, {
        jc1 <- coda.samples(model=jm,
                            variable.names=c(jp, "deviance"),
                            n.iter=1500)
        return(as.mcmc(jc1))
    })
}) ## 13766s


mc1 <- as.mcmc.list(jc1.out)

plot(mc1, ask=TRUE)

summary(mc1[,6:10])


if(!dir.exists("testS")) dir.create("testS")

save(jc1.out, file="testS/jc1.gzip")





## Draw posterior samples
system.time({
    jc1za.out <- clusterEvalQ(cl1, {
        jc1 <- coda.samples(model=jm,
                            variable.names=c(jp, "z", "a"),
                            n.iter=5000)
        return(as.mcmc(jc1))
    })
}) ## 13766s


varnames(as.mcmc.list(jc1za.out))


save(jc1za.out, file="testS/jc1za.gzip")



stopCluster(cl1)



## Model with hole-free state-space

library(parallel)

### Put the data in a list
jd2 <- list(y=y, y.zero=matrix(0, M, nYears),
            z=zdata.aug,
            nBears=nBears, M=M,
            nTrapsYr=as.integer(nTrapsYr),
            nWeeks=nWeeks, nYears=nYears, nFuture=0,
            oper.traps=oper.traps,
            prevcap=prevcap.f+1,
            S=S2, delta=delta,
            ones=matrix(1, M, nYears),
            x=traps, xlim=xlim, ylim=ylim)
jd2$z=cbind(zdata.aug, matrix(NA, M, jd2$nFuture))

str(jd2)





(nCores <- max(1, min(detectCores()-1, 4)))

cl2 <- makeCluster(nCores)

clusterExport(cl2, c("jd1", "jd2", "ji", "jp",
                     "traps.in", "first.det", "last.det"))

system.time({
    jm2.out <- clusterEvalQ(cl2, {
        library(rjags)
        load.module("dic")
        load.module("glm")
        jm2 <- jags.model(file="dT-s0-rPDD2-phi0-pYB-sig0.jag",
                          data=jd2, inits=ji,
                          n.chains=1, n.adapt=500)
        return(jm2)
    })
})


## save(jm.out, file="jm21.gzip")


system.time({
    jc2.1.out <- clusterEvalQ(cl2, {
        jc2.1 <- coda.samples(model=jm2,
                              variable.names=c(jp, "deviance"),
                              n.iter=1500)
        return(as.mcmc(jc2.1))
    })
})


mc2 <- as.mcmc.list(jc2.1.out)

plot(mc2, ask=TRUE)

varnames(mc2)

summary(mc2[,6:10])


if(!dir.exists("testS")) dir.create("testS")

save(jc2.1.out, file="testS/jc2.gzip")




## Draw posterior samples
system.time({
    jc2za.out <- clusterEvalQ(cl2, {
        jc2 <- coda.samples(model=jm2,
                            variable.names=c(jp, "z", "a"),
                            n.iter=5000)
        return(as.mcmc(jc2))
    })
}) ## 13766s


varnames(as.mcmc.list(jc2za.out))


save(jc2za.out, file="testS/jc2za.gzip")



ls()












### Compare abundance estimates

library(coda)

load("testS/jc1za.gzip")
load("testS/jc2za.gzip")


mc1za.N <- as.mcmc.list(jc1za.out)[,6:10]
mc2za.N <- as.mcmc.list(jc2za.out)[,6:10]

Nq1 <- summary(mc1za.N)$quantile[,c(1,3,5)]
Nq2 <- summary(mc2za.N)$quantile[,c(1,3,5)]

Nq1[,2]-Nq2[,2]

pt.pch <- (1:6)[-3]

pdf("../figs/N_SS1vSS2.pdf", width=6, height=6)
plot(Nq1[,2], Nq2[,2], type="n", ##pch=1:5,
     xlab="Abundance (Original state-space)",
     ylab="Abundance (Filled state-space)",
     xlim=c(50, 250), ylim=c(50, 250))
abline(a=0,b=1,lty=2)
segments(Nq2[,2], Nq1[,1], Nq2[,2], Nq1[,3])
segments(Nq2[,1], Nq1[,2], Nq2[,3], Nq1[,2])
points(Nq2[,2], Nq1[,2], pch=pt.pch, cex=1)
legend(50, 250, paste("Year", 2012:2016), pch=pt.pch, pt.cex=1)
dev.off()
system("open ../figs/N_SS1vSS2.pdf")






### ----------------- Forecast with each state-space ----------------



source("forecast_fn.R")


## No harvest (state space 2)
load("testS/jc1za.gzip")

mc1za <- window(as.mcmc.list(jc1za.out), thin=10)

system.time({
    forc1.h0 <- forecast(mc1za, nFuture=50, M=15000,
                         harvest=0, stochastic=TRUE,
                         maxRecruit=1,
                         NbarSims=100,
                         report=100)
})


str(forc2.h0)


system.time({
    forc1.h5 <- forecast(mc1za, nFuture=50, M=15000,
                         harvest=5, stochastic=TRUE,
                         maxRecruit=1,
                         NbarSims=100,
                         report=100)
})


system.time({
    forc1.h10 <- forecast(mc1za, nFuture=50, M=15000,
                          harvest=10, stochastic=TRUE,
                          maxRecruit=1,
                          NbarSims=100,
                          report=100)
})


system.time({
    forc1.h15 <- forecast(mc1za, nFuture=50, M=15000,
                          harvest=15, stochastic=TRUE,
                          maxRecruit=1,
                          NbarSims=100,
                          report=100)
})


system.time({
    forc1.h20 <- forecast(mc1za, nFuture=50, M=15000,
                          harvest=20, stochastic=TRUE,
                          maxRecruit=1,
                          NbarSims=100,
                          report=100)
})


str(forc1.h20)

## Percent of simulations with M too low
mean(forc1.h0$bad)*100
mean(forc1.h5$bad)*100
mean(forc1.h10$bad)*100
mean(forc1.h15$bad)*100
mean(forc1.h20$bad)*100



forc1.h0.Nbar <- forc1.h0$Nbar
forc1.h5.Nbar <- forc1.h5$Nbar
forc1.h10.Nbar <- forc1.h10$Nbar
forc1.h15.Nbar <- forc1.h15$Nbar
forc1.h20.Nbar <- forc1.h20$Nbar

save(forc1.h0.Nbar, forc1.h5.Nbar, forc1.h10.Nbar,
     forc1.h15.Nbar, forc1.h20.Nbar,
     file="testS/forecasts_ss1.gzip")


## No harvest (state space 2)
load("testS/jc2za.gzip")

mc2za <- window(as.mcmc.list(jc2za.out), thin=10)

system.time({
    forc2.h0 <- forecast(mc2za, nFuture=50, M=15000,
                         harvest=0, stochastic=TRUE,
                         maxRecruit=1,
                         NbarSims=100,
                         report=100)
})


str(forc2.h0)


system.time({
    forc2.h5 <- forecast(mc2za, nFuture=50, M=15000,
                         harvest=5, stochastic=TRUE,
                         maxRecruit=1,
                         NbarSims=100,
                         report=100)
})


str(forc2.h5)




system.time({
    forc2.h10 <- forecast(mc2za, nFuture=50, M=15000,
                          harvest=10, stochastic=TRUE,
                          maxRecruit=1,
                          NbarSims=100,
                          report=100)
})

str(forc2.h10)


system.time({
    forc2.h15 <- forecast(mc2za, nFuture=50, M=15000,
                          harvest=15, stochastic=TRUE,
                          maxRecruit=1,
                          NbarSims=100,
                          report=100)
})

str(forc2.h15)




system.time({
    forc2.h20 <- forecast(mc2za, nFuture=50, M=15000,
                         harvest=20, stochastic=TRUE,
                         maxRecruit=1,
                         NbarSims=100,
                         report=100)
})


str(forc2.h20)


## Percent of simulations with M too low
mean(forc2.h0$bad)*100
mean(forc2.h5$bad)*100
mean(forc2.h10$bad)*100
mean(forc2.h15$bad)*100
mean(forc2.h20$bad)*100


forc2.h0.Nbar <- forc2.h0$Nbar
forc2.h5.Nbar <- forc2.h5$Nbar
forc2.h10.Nbar <- forc2.h10$Nbar
forc2.h15.Nbar <- forc2.h15$Nbar
forc2.h20.Nbar <- forc2.h20$Nbar

save(forc2.h0.Nbar, forc2.h5.Nbar, forc2.h10.Nbar,
     forc2.h15.Nbar, forc2.h20.Nbar,
     file="testS/forecasts_ss2.gzip")








## Compare extinction risk in two state-spaces

load("testS/forecasts_ss1.gzip")
load("testS/forecasts_ss2.gzip")

EX1.h0 <- apply(forc1.h0.Nbar==0, c(1, 3), mean)
EX1.h5 <- apply(forc1.h5.Nbar==0, c(1, 3), mean)
EX1.h10 <- apply(forc1.h10.Nbar==0, c(1, 3), mean)
EX1.h15 <- apply(forc1.h10.Nbar==0, c(1, 3), mean)
EX1.h20 <- apply(forc1.h20.Nbar==0, c(1, 3), mean)

EX2.h0 <- apply(forc2.h0.Nbar==0, c(1, 3), mean)
EX2.h5 <- apply(forc2.h5.Nbar==0, c(1, 3), mean)
EX2.h10 <- apply(forc2.h10.Nbar==0, c(1, 3), mean)
EX2.h15 <- apply(forc2.h10.Nbar==0, c(1, 3), mean)
EX2.h20 <- apply(forc2.h20.Nbar==0, c(1, 3), mean)


EX1.h0.mean <- rowMeans(EX1.h0)
EX1.h0.low <- apply(EX1.h0, 1, quantile, prob=0.025)
EX1.h0.upp <- apply(EX1.h0, 1, quantile, prob=0.975)
EX1.h5.mean <- rowMeans(EX1.h5)
EX1.h5.low <- apply(EX1.h5, 1, quantile, prob=0.025)
EX1.h5.upp <- apply(EX1.h5, 1, quantile, prob=0.975)
EX1.h10.mean <- rowMeans(EX1.h10)
EX1.h10.low <- apply(EX1.h10, 1, quantile, prob=0.025)
EX1.h10.upp <- apply(EX1.h10, 1, quantile, prob=0.975)
EX1.h15.mean <- rowMeans(EX1.h15)
EX1.h15.low <- apply(EX1.h15, 1, quantile, prob=0.025)
EX1.h15.upp <- apply(EX1.h15, 1, quantile, prob=0.975)
EX1.h20.mean <- rowMeans(EX1.h20)
EX1.h20.low <- apply(EX1.h20, 1, quantile, prob=0.025)
EX1.h20.upp <- apply(EX1.h20, 1, quantile, prob=0.975)


EX2.h0.mean <- rowMeans(EX2.h0)
EX2.h0.low <- apply(EX2.h0, 1, quantile, prob=0.025)
EX2.h0.upp <- apply(EX2.h0, 1, quantile, prob=0.975)
EX2.h5.mean <- rowMeans(EX2.h5)
EX2.h5.low <- apply(EX2.h5, 1, quantile, prob=0.025)
EX2.h5.upp <- apply(EX2.h5, 1, quantile, prob=0.975)
EX2.h10.mean <- rowMeans(EX2.h10)
EX2.h10.low <- apply(EX2.h10, 1, quantile, prob=0.025)
EX2.h10.upp <- apply(EX2.h10, 1, quantile, prob=0.975)
EX2.h15.mean <- rowMeans(EX2.h15)
EX2.h15.low <- apply(EX2.h15, 1, quantile, prob=0.025)
EX2.h15.upp <- apply(EX2.h15, 1, quantile, prob=0.975)
EX2.h20.mean <- rowMeans(EX2.h20)
EX2.h20.low <- apply(EX2.h20, 1, quantile, prob=0.025)
EX2.h20.upp <- apply(EX2.h20, 1, quantile, prob=0.975)


years.f <- seq(2016, length.out=51)

pdf("../figs/EX_SS1vSS2.pdf", width=7, height=7)
par(mfrow=c(3,2), mai=c(0.6,0.6,0.3,0.1), omi=c(0,0,0,0))
plot(years.f, EX1.h0.mean, type="l", col=1, ylim=c(0, 0.4),
     main="No increased harvest",
     xlab="Year", ylab="Extinction risk", lwd=2)
lines(years.f, EX2.h0.mean, type="l", col=3)
## lines(EX3.h0.mean, type="l", col=3, ylim=c(0, 0.5))
legend(2016, 0.4, c("Original state-space", "Filled state-space"),
       col=c(1,3), lty=1, lwd=c(2,1))
plot(years.f, EX1.h5.mean, type="l", col=1, ylim=c(0, 0.4),
     main="Harvest increased by 5 females/year",
     xlab="Year", ylab="Extinction risk", lwd=2)
lines(years.f, EX2.h5.mean, type="l", col=3)
## lines(EX3.h5.mean, type="l", col=3, ylim=c(0, 0.5))
plot(years.f, EX1.h10.mean, type="l", col=1, ylim=c(0, 0.4),
     main="Harvest increased by 10 females/year",
     xlab="Year", ylab="Extinction risk", lwd=2)
lines(years.f, EX2.h10.mean, type="l", col=3)
## lines(EX3.h10.mean, type="l", col=3, ylim=c(0, 0.5))
plot(years.f, EX1.h15.mean, type="l", col=1, ylim=c(0, 0.4),
     main="Harvest increased by 15 females/year",
     xlab="Year", ylab="Extinction risk", lwd=2)
lines(years.f, EX2.h15.mean, type="l", col=3)
## lines(EX3.h15.mean, type="l", col=3, ylim=c(0, 0.5))
plot(years.f, EX1.h20.mean, type="l", col=1, ylim=c(0, 0.4),
     main="Harvest increased by 20 females/year",
     xlab="Year", ylab="Extinction risk", lwd=2)
lines(years.f, EX2.h20.mean, type="l", col=3)
## lines(EX3.h20.mean, type="l", col=3, ylim=c(0, 0.5)
dev.off()
system("open ../figs/EX_SS1vSS2.pdf")


plot(EX2.h20.mean, EX1.h20.mean)
abline(a=0,b=1,lty=2)

plot(EX2.h10.mean, EX1.h10.mean)
abline(a=0,b=1,lty=2)







