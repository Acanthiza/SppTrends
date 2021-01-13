
### 	Set working directory
setwd("C:/Workspace_NW/LRG/R")

###	Calculate Extinction Probability

#	Import and setup data
rm(list=ls())
require(RODBC)
con <- odbcConnectAccess("C:/Workspace_NW/LRG/rel/LRG_rel.mdb")
dat <- sqlFetch(con,"q01650")
odbcCloseAll()
rm(con)

head(dat)

names(dat) <- c("ls","type","ttype","com","n","run","source","stat")

#	After importing data, change the ls== in the next call

dat.sub <- subset(dat,ls=="NORTH HILLY FLINDERS")

head(dat.sub)

library(reshape)
dat.cast <- cast(dat.sub, ls + com ~ source, mean, fill=0.5, value="stat")

#	Remove NA column, if it exists (depends on whether any species are not indicated as declining by any method)
dat.cast <- if (colnames(dat.cast)[max(length(colnames(dat.cast)))] == "NA") dat.cast[,-max(length(colnames(dat.cast)))] else dat.cast

summary(dat.cast)

ls <- dat.cast$ls[1]

dat.pca <- princomp(dat.cast[,-c(1:2)])
#, scale = TRUE

plot(dat.pca)
loadings(dat.pca)

lab = dat.cast$com
lab.2 = colnames(dat.cast[,-c(1:2)])

windowsFonts(A=windowsFont("Century Gothic"))
par(mfrow=c(1,1),family="A",mar=c(4,4,5,2)+0.1)
biplot(dat.pca$scores[,c(1,2)],dat.pca$loadings[,c(1,2)], xlabs = lab, ylabs = lab.2,cex=0.5)
#biplot(dat.pca$scores[,c(1,3)],dat.pca$loadings[,c(1,3)], xlabs = lab, ylabs = lab.2,cex=0.5)
title(dat.cast$ls[1])


#	RUN THROUGH TO HERE - then subjectively select cutoff scores for decline and put into the for statement below

pca.dec <- data.frame(ls=ls,com=dat.cast$com,dat.pca$scores)

for (i in 1:length(row.names(pca.dec))){
	pca.dec$pcadec[i] <- if (pca.dec$Comp.1[i]> 0.4 & pca.dec$Comp.2[i]< 0.14) "dec" else "notdec"								# one y criteria
#	pca.dec$pcadec[i] <- if (pca.dec$Comp.1[i]> 0.2 & pca.dec$Comp.2[i]> -0.15 | pca.dec$com[i] == "Brush Bronzewing") "dec" else "notdec"	# two y criteria
#	pca.dec$pcadec[i] <- if (pca.dec$Comp.1[i]> 0.25 & pca.dec$Comp.3[i]< 0) "dec" else "notdec"								# one y criteria on Comp.3
#	pca.dec$pcadec[i] <- if (pca.dec$Comp.1[i]> 0.35 | pca.dec$Comp.2[i] == 0.430728471) "dec" else "notdec"						# random
}

head(pca.dec)

##	Print pca plots
#	Dec vs. not dec
png(file=paste("C:/Workspace_NW/LRG/R/pcadec/",ls,"dec.png",sep=""),width=2000,height=2000,res=200)
windowsFonts(A=windowsFont("Century Gothic"))
par(family="A", mar=c(4,4,5,2)+0.1)
biplot(dat.pca$scores[,c(1,2)],dat.pca$loadings[,c(1,2)], xlabs = pca.dec$pcadec, ylabs = lab.2,cex=0.7)
title(dat.cast$ls[1])
dev.off()

#	Spp
png(file=paste("C:/Workspace_NW/LRG/R/pcadec/",ls,"spp.png",sep=""),width=2000,height=2000,res=200)
windowsFonts(A=windowsFont("Century Gothic"))
par(family="A", mar=c(4,4,5,2)+0.1)
biplot(dat.pca$scores[,c(1,2)],dat.pca$loadings[,c(1,2)], xlabs = lab, ylabs = lab.2,cex=0.7)
title(dat.cast$ls[1])
dev.off()


#	Once happy with result, append results to .csv

# export results
write.table(pca.dec,file="pcadec.csv",row.names=F,append=T,sep=",",col.names=F)






