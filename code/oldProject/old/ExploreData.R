

##Import and setup data
rm(list=ls())
require(RODBC)
con <- odbcConnectAccess("LRG/rel/LRG_rel.mdb")
dat <- sqlFetch(con,"q01120")
odbcCloseAll()
rm(con)

head(dat)

names(dat) <- c("poly","type","spp","year","hex","n")

dat$decade <- round(dat$year/10,0)*10
dat$decade <- as.factor(dat$decade)

library(reshape)

dat.cast <- cast(dat,poly + type + spp ~ decade, sum, value = "n")
head(dat.cast)

sub.dat <- subset(dat.cast, poly=="CENTRAL HILLS")
head(sub.dat)

library(vegan)

dat.mds <- metaMDS(sub.dat[,-c(1:3)])

lab <- dat.cast$decade

windowsFonts(A=windowsFont("Century Gothic"))
par(mfrow=c(1,1),family="A")
plot(dat.mds$points[,1], dat.mds$points[,2],xlab='Dim1',ylab='Dim2',type="n")
text(dat.mds$points[,1], dat.mds$points[,2],labels=lab,cex=0.5)
text(dat.mds$species[,1], jitter(dat.mds$species[,2]), labels=colnames(dat.cast), cex=0.8, col="red")
#text(dat.mds$points[,1], dat.mds$points[,2],cex=0.5)# to see mds with ID number instead of group number
title(paste("Stress = ",round(dat.mds$stress,2)," for ",dat.mds$ndim,"-d solution." ,sep=""),cex=1)
ordihull(dat.mds, cut)

