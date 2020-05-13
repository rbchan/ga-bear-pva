## Fit open population SCR models to central GA black bear data
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
y.nonsp.nonrobust <- apply(y.nonsp.robust>0,
                           c(1,3), any)*1L
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
## Models are slow. About 450 iter/hr

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
    out <- list(z=zi, s=si, s.rec=si.rec, p0=p0i,
                sigma=runif(1, 3000, 4000), behave=runif(1, 0.5, 1),
                mu.lphi=0, sig.lphi=.01,
                psi=c(runif(1), rep(NA, nYears+nFuture-1)),
                igamma0=runif(1, 0.2, 0.3), gamma1=runif(1, 0, 0.001))
    out$".RNG.name" <- "base::Mersenne-Twister" 
    out$".RNG.seed" <- sample.int(1e6, 1)
    return(out)
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
jp <- c("p0", "sigma", "behave", 
        "igamma0", "gamma0", "gamma1", "phi", "gamma", "EN", 
        "N", "Ntotal", "survivors", "recruits", "lambda")


## Models in parallel with WAIC computation

library(parallel)

(nCores <- max(1, min(detectCores()-1, 4)))

cl1 <- makeCluster(nCores)

clusterExport(cl1, c("jd", "ji", "jp",
                     "traps.in", "first.det", "last.det"))

## Compile and adapt. Comment out next line if you don't want reproducible results
clusterSetRNGStream(cl1, iseed=389761)

system.time({
    jm.out <- clusterEvalQ(cl1, {
        library(rjags)
        load.module("dic")
        load.module("glm")
        jm <- jags.model(file="dT-s0-rPDD2-phi0-pYB-sig0.jag",
                         data=jd, inits=ji,
                         n.chains=1, n.adapt=500)
        return(jm)
    })
}) ## Approx 1.5 hr


save(jm.out, file="jm.gzip")


## Draw posterior samples, but don't monitor latent variables
system.time({
    jc1.out <- clusterEvalQ(cl1, {
        jc1 <- coda.samples(model=jm,
                            variable.names=c(jp, "deviance"),
                            n.iter=5000)
        return(as.mcmc(jc1))
    })
}) ## 450 iterations/hr


save(jc1.out, file="jc1.gzip")


## mc1 <- as.mcmc.list(jc1.out)
## plot(mc1, ask=TRUE)

## Draw some more samples, including the latent varialbes needed
## for forecasting
jc2.out <- clusterEvalQ(cl1, {
    jc2 <- coda.samples(model=jm,
                        variable.names=c(jp, "deviance", "z", "a","s"),
                        n.iter=5000, n.thin=1)
    return(as.mcmc(jc2))
})



save(jc2.out, file="jc2.gzip")

rm(jc1.out, jc2.out)
gc()


## WAIC
js3.out <- clusterEvalQ(cl1, {
    js3 <- jags.samples(model=jm,
                        variable.names=c("deviance","WAIC"),
                        n.iter=5000, n.thin=1, type="mean")
    return(js3)
})


save(js3.out, file="js3.gzip")



## These additional samples aren't really necessary because the Markov
## chains usually converge very rapidly
jc4.out <- clusterEvalQ(cl1, {
    jc4 <- coda.samples(model=jm,
                        variable.names=c(jp, "deviance", "z", "a","s"),
                        n.iter=5000, n.thin=1)
    return(as.mcmc(jc4))
})



save(jc4.out, file="jc4.gzip")




## WAIC


js5.out <- clusterEvalQ(cl1, {
    js5 <- jags.samples(model=jm,
                        variable.names=c("deviance","WAIC"),
                        n.iter=1000, n.thin=1, type="mean")
    return(js5)
})


save(js5.out, file="js5.gzip")




jm.out <- clusterEvalQ(cl1, {
    return(jm)
})


save(jm.out, file="jm.gzip")



state.out <- clusterEvalQ(cl1, {
    return(jm$state())
})


save(state.out, file="jm-state.gzip")






















ls()

rm(list=ls())
gc()


library(coda)


## load("jc1.gzip")
## jc1 <- as.mcmc.list(jc1.out)

## load("jc2.gzip")
## jc2 <- as.mcmc.list(jc2.out)

load("jc4.gzip")
jc4 <- as.mcmc.list(jc4.out)



vn <- varnames(jc1)

vn.sub <- c(vn[-c(grep("p0", vn), grep("z\\[", vn),
                  grep("s\\[", vn), grep("a\\[", vn))],
            "p0[1,1,1]", "p0[1,1,2]", "p0[1,2,1]", "p0[1,2,2]",
            "p0[1,3,1]", "p0[1,3,2]", "p0[1,4,1]", "p0[1,4,2]",
            "p0[1,5,1]", "p0[1,5,2]",
            paste("gamma[", 1:4, "]", sep=""),
            paste("survivors[", 1:4, "]", sep=""),
            paste("recruits[", 1:4, "]", sep=""),
            paste("lambda[", 1:4, "]", sep=""), "sigma")

##vn21sub <- vn21sub[-grep("deviance", vn21sub)]

mc.sub <- as.data.frame(as.matrix(jc4[,vn.sub]))
summary(mc.sub)

## Area of the state-space
## sum(jdata15$S)*jdata15$delta^2/1e6
SSarea <- 622.5984 ## sq-km



## Expected density (bears per 100 sq-km)
mc.sub$"ED[1]" <- mc.sub$"EN[1]"/SSarea*100
mc.sub$"ED[2]" <- mc.sub$"EN[2]"/SSarea*100
mc.sub$"ED[3]" <- mc.sub$"EN[3]"/SSarea*100
mc.sub$"ED[4]" <- mc.sub$"EN[4]"/SSarea*100
mc.sub$"ED[5]" <- mc.sub$"EN[5]"/SSarea*100

## ## Baseline cap probs
## mc.sub$p0pre <- mc.sub$"p0[1,1,1]"
## mc.sub$p0post <- mc.sub$"p0[1,1,1]" / mc.sub$behave

summary(mc.sub)


mc.sub <- data.matrix(mc.sub)



sstats.0 <- cbind(Mean=colMeans(mc.sub),
                  SD=apply(mc.sub, 2, sd),
                  LowerCI=apply(mc.sub, 2, quantile, prob=0.025),
                  Median=apply(mc.sub, 2, quantile, prob=0.5),
                  UpperCI=apply(mc.sub, 2, quantile, prob=0.975))

vn.out <- c(paste("N[", 1:5, "]", sep=""),
            paste("ED[", 1:5, "]", sep=""),
            "Ntotal",
            "phi",
            paste("survivors[", 1:4, "]", sep=""),
            "igamma0", "gamma1",
            paste("gamma[", 1:4, "]", sep=""),
            paste("recruits[", 1:4, "]", sep=""),
            paste("lambda[", 1:4, "]", sep=""),
            ##             "p0pre", "p0post",
            grep("p0", rownames(sstats.0), value=TRUE),
            "sigma")




sstats <- sstats.0[vn.out,]


sstats

write.table(format(sstats, digits=1, sci=FALSE),
            quote=FALSE, sep="\t")

write.table(format(sstats, digits=1, sci=FALSE),
            file="sstats21-2.txt",
            quote=FALSE, sep="\t")


write.table(round(sstats, digits=2),
            quote=FALSE, sep="\t")

write.table(round(sstats, digits=2),
            file="sstats21-2.txt",
            quote=FALSE, sep="\t")





## Density-dependent recruitment

gamma01 <- as.matrix(jc4[,c("igamma0","gamma1")])[seq(1, 20000, by=10),]

str(gamma01)

Nx <- seq(0, 200)

DDpost <- apply(gamma01, 1, function(x) x[1]*exp(-x[2]*Nx^2))

DDpost.mean <- apply(DDpost, 1, mean)
DDpost.low <- apply(DDpost, 1, quantile, prob=0.025)
DDpost.upp <- apply(DDpost, 1, quantile, prob=0.975)


png("../figs/fig-3_dd-recruitment.png", width=6, height=5,
    units="in", res=500)
##pdf("fig/fig-3_dd-recruitment.pdf", width=6, height=5)
par(mai=c(0.8, 0.8, 0.1, 0.1))
plot(sstats.0[paste("N[", 1:4, "]", sep=""),"Median"],
     sstats.0[paste("gamma[", 1:4, "]", sep=""),"Median"],
     xlim=c(0, 200), xaxt="n",
     ylim=c(0, 1), xaxs="i", yaxs="i",
     frame=FALSE, 
     xlab=expression(paste("Density (females / 100 ", km^2, ")")),
##     xlab="Female abundance", 
     ylab=expression(paste(italic("Per capita"), " recruitment")))
##matlines(DDpost, type="l", col=gray(0.8))
points(sstats.0[paste("N[", 1:4, "]", sep=""),"Median"],
       sstats.0[paste("gamma[", 1:4, "]", sep=""),"Median"], pch=16)
arrows(sstats.0[paste("N[", 1:4, "]", sep=""),"Median"],
       sstats.0[paste("gamma[", 1:4, "]", sep=""),"LowerCI"],
       sstats.0[paste("N[", 1:4, "]", sep=""),"Median"],
       sstats.0[paste("gamma[", 1:4, "]", sep=""),"UpperCI"],
       angle=90, code=3, length=0.05)
lines(Nx, DDpost.mean)
lines(Nx, DDpost.low, lty=2)
lines(Nx, DDpost.upp, lty=2)
Dx <- seq(0, 40, 10)
##axis(1, at=seq(0, 200, 50), labels=round(seq(0, 200, 50)/SSarea*100, 1))
axis(1, at=Dx*SSarea/100, labels=Dx)
dev.off()
system("open ../figs/fig-3_dd-recruitment.png")
##system("open ../figs/fig-3_dd-recruitment.pdf")




## Annual abundance and density estimates
png("../figs/fig-5_abundance.png", width=8, height=6, units="in", res=400)
par(mai=c(0.9, 0.9, 0.9, 0.9))
yrs <- 2012:2016
plot(yrs, sstats.0[paste("N[", 1:5, "]", sep=""),"Median"],
     type="b", ylim=c(70, 210), pch=16,
     xlab="Year", ylab="Female abundance", cex.lab=1.5)
axis(4, at=seq(80, 200, 20),
     label=round(seq(80, 200, 20)/SSarea*100, 1))
mtext(expression(paste("Density (females / 100 ", km^2, ")")),
      side=4, line=-1.5,
      outer=TRUE, cex=1.5)
arrows(yrs,
       sstats.0[paste("N[", 1:5, "]", sep=""),"LowerCI"],
       yrs,
       sstats.0[paste("N[", 1:5, "]", sep=""),"UpperCI"],
       angle=90, code=3, length=0.05)
dev.off()
system("open ../figs/fig-5_abundance.png")







