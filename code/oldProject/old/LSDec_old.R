
### 	Set working directory
setwd("C:/Workspace_NW/LRG/R")

##Import and setup data
rm(list=ls())
require(RODBC)
con <- odbcConnectAccess("C:/Workspace_NW/LRG/rel/LRG_rel_old.mdb")
dat <- sqlFetch(con,"q01700")
odbcCloseAll()
rm(con)

head(dat)

names(dat) <- c("poly","type","ttype","dec","notdec","pdec","rem","area")

dat[is.na(dat)] <- 0

dat$spprich <- dat$dec + dat$notdec

for (i in 1:length(dat$rem)){
	if (dat$rem[i] < 0.1) {dat$MandH[i] <- "red"}
	if (dat$rem[i] > 0.1 & dat$rem[i] < 0.4)  {dat$MandH[i] <- "orange"}
	if (dat$rem[i] > 0.4 & dat$rem[i] < 0.9)  {dat$MandH[i] <- "yellow"}
	if (dat$rem[i] > 0.9) {dat$MandH[i] <- "green"}
	}	

###	Linear model of prop decline vs. spprichness (more species, more likely to have more declining?)
dec.mod <- lm(pdec ~ area, data=dat)


###	Create graph
#	Set up a file that is written to later in the script
png(file=paste("C:/Workspace_NW/LRG/R/lsdec/lsdec.png",sep=""),width=3000,height=2000,res=200)

#	Set up graphical parameters
windowsFonts(A=windowsFont("Century Gothic"))
par(mfrow=c(1,1),family="A",mar=c(5,4,2,2)+0.1)

#	Create plot
plot(dat$pdec~dat$area,xlab="Landscape Area (hectares)",
	pch = 20,
	ylab="Proportion of species declining", font.lab=2,
	sub=paste("Red fitted line is linear model (p=",
	round(coef(dec.mod)[2],3),
	") +/- s.e. Circle size indicates landscape remnancy. Circle colours represent Hobbs (2004) landscape classification (see key on chart).",
	sep=""))

#	Add regression line with errors
x.val <- data.frame(area=seq(min(dat$area),max(dat$area),1000))
pred <- predict(dec.mod,x.val,se.fit=T)
seup <- pred$fit+pred$se.fit
sedown <- pred$fit-pred$se.fit
lines(x.val$area,pred$fit,col="red",lwd=2)
lines(x.val$area,seup,col="red",lty=2)
lines(x.val$area,sedown,col="red",lty=2)

#	Add data points (including circles)
symbols(dat$area,dat$pdec,dat$rem,add=T,fg=(dat$MandH), inches = 1)
points(dat$area,dat$pdec,pch=16,col=(dat$MandH))

#	Add labels
require(maptools)
placement <- pointLabel(dat$area, dat$pdec, labels = dat$poly,
	cex = 0.8, doPlot=F)
text(placement, labels=dat$poly, cex=0.8, offset=0)

#	Add legend
legend(0.9*max(dat$area),max(dat$pdec),
	c("relictual","fragmented","variegated","intact"),
	pch=16, col=c("red","orange","yellow","green"),
	bg="grey")

#	Write the plot to file
dev.off()



