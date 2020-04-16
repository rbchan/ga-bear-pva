## Load the data
load("../data/ga_bear_data_females.gzip")
ls()

## Dimensions
bear.ids <- rownames(ch4d.f)
dim.ch4d <- dim(ch4d.f)
nBears <- dim.ch4d[1]
nWeeks <- dim.ch4d[2]
nTraps <- dim.ch4d[3]
nYears <- dim.ch4d[4]

all(rownames(traps)==dimnames(ch4d.f)[[3]]) ## Must be TRUE


## Non-spatial, robust design format
y.nonsp.robust <- apply(ch4d.f>0, c(1,3,4), any)*1L
str(y.nonsp.robust)

## Non-spatial, non-robust format
y.nonsp.nonrobust <- apply(y.nonsp.robust>0, c(1,3), any)*1L
str(y.nonsp.nonrobust)

## Spatial, non-robust format
y.sp.nonrobust <- apply(ch4d.f>0, c(1,2,4), any)*1L
str(y.sp.nonrobust)

## First and last detections
first.det <- apply(y.nonsp.nonrobust>0, 1, function(x) min(which(x)))
last.det <- apply(y.nonsp.nonrobust>0, 1, function(x) max(which(x)))
known.alive <- last.det-first.det+1
table(known.alive)



## Recap stats

table(rowSums(ch4d.f)) ## Total recap frequencies

table(apply(apply(y.sp.nonrobust>0, c(1,2), any), 1, sum)) ## Spatial recaps

apply(apply(y.sp.nonrobust>0, c(1,3), sum), 2, table) ## Spatial recaps by year


sum(y.nonsp.nonrobust[,1])
sum(y.nonsp.nonrobust[,2]*(1-y.nonsp.nonrobust[,1])) # nObservedRecruits
sum(y.nonsp.nonrobust[,3]*(1-y.nonsp.nonrobust[,2])) #
sum(y.nonsp.nonrobust[,4]*(1-y.nonsp.nonrobust[,3])) #


colSums(y.nonsp.nonrobust)

image(t(y.nonsp.nonrobust))



## Assess behavioral response
firstWeekYr1 <- apply(y.nonsp.robust[first.det==1,,1]>0, 1, which.max)
firstWeekYr2 <- apply(y.nonsp.robust[first.det==2,,2]>0, 1, which.max)
firstWeekYr3 <- apply(y.nonsp.robust[first.det==3,,3]>0, 1, which.max)
firstWeekYr4 <- apply(y.nonsp.robust[first.det==4,,4]>0, 1, which.max)
firstWeekYr5 <- apply(y.nonsp.robust[first.det==5,,5]>0, 1, which.max)

table(firstWeekYr1)#/sum(firstWeekYr1))
table(firstWeekYr2)#/sum(firstWeekYr2)
table(firstWeekYr3)#/sum(firstWeekYr3)
table(firstWeekYr4)#/sum(firstWeekYr4)
table(firstWeekYr5)#/sum(firstWeekYr5)


plot(table(firstWeekYr1)/sum(firstWeekYr1))
plot(table(firstWeekYr2)/sum(firstWeekYr2))
plot(table(firstWeekYr3)/sum(firstWeekYr3))
plot(table(firstWeekYr4)/sum(firstWeekYr4))
plot(table(firstWeekYr5)/sum(firstWeekYr5))






## Assess evidence of dispersal

## Average locations for each bear in each year
## Ignore within-year recaps, which are strongly influenced by bait
avloc <- array(NA_real_, c(nBears, nYears, 2))
for(i in 1:nBears) {
    for(t in 1:nYears) {
        y.it <- y.sp.nonrobust[i,,t]
        if(all(y.it==0))
            next
        traps.it <- traps[which(y.it>0),,drop=FALSE]
        avloc[i,t,] <- colMeans(traps.it)
    }
}

avdist <- list() ##vector(mode="list", length=nBears)
counter <- 1
for(i in 1:nBears) {
    if(sum(!is.na(avloc[i,,1]))<2) {
        next
    }
    avdist[[counter]] <- dist(avloc[i,,])
    counter <- counter+1
}

avdist.means <- sapply(avdist, mean, na.rm=TRUE)
avdist.maxs <- sapply(avdist, max, na.rm=TRUE)
mean(avdist.means)
range(avdist.means)
sd(avdist.means)

mean(avdist.maxs)
max(avdist.maxs)













## Map bear detections

if(!dir.exists("../figs"))
    dir.create("../figs")
if(!dir.exists("../figs/detmaps"))
    dir.create("../figs/detmaps")

clrs <- col2rgb(c("black", "blue", "green", "red", "tan"), TRUE)
clrs[4,] <- 125
clrs <- rgb(red=clrs[1,], green=clrs[2,], blue=clrs[3,],
            alpha=clrs[4,], maxColorValue=255)

for(i in 1:nBears) {
    cat("Doing bear", i, "of", nBears, "\n")
    file.i <- paste("../figs/detmaps/bear", bear.ids[i],
                    ".png", sep="")
    png(file.i, width=5, height=6, units="in", res=400)
    par(mai=c(0.5, 0.1, 0.5, 0.1))
    plot(traps/1000, pch=3, asp=1, axes=FALSE, frame=FALSE,
         xlab="", ylab="", main=bear.ids[i])
    legend("topleft", legend=2012:2016, pch=21:25, pt.bg=clrs)
    axis(1, at=c(248, 253, 258), label=c("0", "", "10 km"))
    for(t in 1:nYears) {
        if(y.nonsp.nonrobust[i,t]<1)
            next
        detlocs.it <- y.sp.nonrobust[i,,t]==1
        points(traps[detlocs.it,,drop=FALSE]/1000, pch=20+t,
               bg=clrs[t], cex=1.5-(t-1)*0.2)
    }
    dev.off()
}



## Trap dets in each year

## 2012
png("../figs/dets2012.png", width=5, height=6, units="in", res=400)
par(mai=c(0.5, 0.1, 0.5, 0.1))
plot(traps/1000, pch=3, asp=1, axes=FALSE, frame=FALSE,
         xlab="", ylab="", main="2012")
points(traps/1000, cex=colSums(y.sp.nonrobust[,,1]),
       pch=16, col=rgb(0,0,1,0.5))
legend("topleft", legend=1:4, pch=16, pt.cex=1:4,
       col=rgb(0,0,1,0.5), title="Detections", y.intersp=1.5)
axis(1, at=c(248, 253, 258), label=c("0", "", "10 km"))
dev.off()
system("open ../figs/dets2012.png")


## 2013
png("../figs/dets2013.png", width=5, height=6, units="in", res=400)
par(mai=c(0.5, 0.1, 0.5, 0.1))
plot(traps/1000, pch=3, asp=1, axes=FALSE, frame=FALSE,
         xlab="", ylab="", main="2013")
points(traps/1000, cex=colSums(y.sp.nonrobust[,,2]),
       pch=16, col=rgb(0,0,1,0.5))
legend("topleft", legend=1:4, pch=16, pt.cex=1:4,
       col=rgb(0,0,1,0.5), title="Detections", y.intersp=1.5)
axis(1, at=c(248, 253, 258), label=c("0", "", "10 km"))
dev.off()
system("open ../figs/dets2013.png")


## 2014
png("../figs/dets2014.png", width=5, height=6, units="in", res=400)
par(mai=c(0.5, 0.1, 0.5, 0.1))
plot(traps/1000, pch=3, asp=1, axes=FALSE, frame=FALSE,
         xlab="", ylab="", main="2014")
points(traps/1000, cex=colSums(y.sp.nonrobust[,,3]),
       pch=16, col=rgb(0,0,1,0.5))
legend("topleft", legend=1:4, pch=16, pt.cex=1:4,
       col=rgb(0,0,1,0.5), title="Detections", y.intersp=1.5)
axis(1, at=c(248, 253, 258), label=c("0", "", "10 km"))
dev.off()
system("open ../figs/dets2014.png")


## 2015
png("../figs/dets2015.png", width=5, height=6, units="in", res=400)
par(mai=c(0.5, 0.1, 0.5, 0.1))
plot(traps/1000, pch=3, asp=1, axes=FALSE, frame=FALSE,
         xlab="", ylab="", main="2015")
points(traps/1000, cex=colSums(y.sp.nonrobust[,,4]),
       pch=16, col=rgb(0,0,1,0.5))
legend("topleft", legend=1:4, pch=16, pt.cex=1:4,
       col=rgb(0,0,1,0.5), title="Detections", y.intersp=1.5)
axis(1, at=c(248, 253, 258), label=c("0", "", "10 km"))
dev.off()
system("open ../figs/dets2015.png")


## 2016
png("../figs/dets2016.png", width=5, height=6, units="in", res=400)
par(mai=c(0.5, 0.1, 0.5, 0.1))
plot(traps/1000, pch=3, asp=1, axes=FALSE, frame=FALSE,
     xlab="", ylab="", main="2016")
points(traps/1000, cex=colSums(y.sp.nonrobust[,,5]),
       pch=16, col=rgb(0,0,1,0.5))
legend("topleft", legend=1:4, pch=16, pt.cex=1:4,
       col=rgb(0,0,1,0.5), title="Detections", y.intersp=1.5)
axis(1, at=c(248, 253, 258), label=c("0", "", "10 km"))
dev.off()
system("open ../figs/dets2016.png")

