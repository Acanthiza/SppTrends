#
#Import and setup data
#

	rm(list=ls())
	require(RODBC)
	con <- odbcConnectAccess("LRG/rel/LRG_rel.mdb")
	dat <- sqlFetch(con,"q01700")
	odbcCloseAll()
	rm(con)

	head(dat)

	names(dat) <- c("poly","plantspprich","type","ttype","dec","notdec","pdec","rem","area")

	dat[is.na(dat)] <- 0

	dat$spprich <- dat$dec + dat$notdec

	for (i in 1:length(dat$rem)){
		if (dat$rem[i] < 0.1) {dat$MandH[i] <- "red"}
		if (dat$rem[i] > 0.1 & dat$rem[i] < 0.4)  {dat$MandH[i] <- "orange"}
		if (dat$rem[i] > 0.4 & dat$rem[i] < 0.9)  {dat$MandH[i] <- "yellow"}
		if (dat$rem[i] > 0.9) {dat$MandH[i] <- "green"}
	}	

#
# Linear model of prop decline vs. spprichness (more species, more likely to have more declining?)
#

	dec.mod <- lm(pdec ~ spprich, data=dat)

#
# Ensure directory for figure outputs exists
#

	dir.create(file.path("LRG","figs","LSDec"))

#
# Create graph
#

	png(file=file.path("LRG","figs","LSDec","lsdec.png"),width=3000,height=2000,res=200)

#
# Set up graphical parameters
#

	windowsFonts(A=windowsFont("Century Gothic"))
	par(mfrow=c(1,1),family="A",mar=c(5,4,2,2)+0.1)

#
# Create plot
#

	plot(dat$pdec~dat$spprich,xlab="Landscape Species Richness"
		, pch = 20
		, ylab="Proportion of species declining", font.lab=2
		, sub=paste("Red fitted line is linear model (p="
		, round(coef(dec.mod)[2],3)
		, ") +/- s.e. Circle size indicates landscape area. Circle colours represent Hobbs (2004) landscape classification (see key on chart).",
		, sep="")
		)

#
# Add data points (including circles)
#

	symbols(dat$spprich,dat$pdec,dat$area,add=T,bg=(dat$MandH),col=(dat$MandH))
	points(dat$spprich,dat$pdec,pch=1,cex=0.5)

#
# Add labels
#

	require(maptools)
	placement <- pointLabel(dat$spprich, dat$pdec, labels = dat$poly
		, cex = 0.8, doPlot=F)

	text(placement, labels=dat$poly, cex=0.8, offset=0)

#
# Add legend
#

	legend(min(dat$spprich),max(dat$pdec),
		c("relictual","fragmented","variegated","intact"),
		pch=19, col=c("red","orange","yellow","green"),
		bg="grey")

#
# Add regression line with errors
#

	x.val <- data.frame(spprich=seq(min(dat$spprich),max(dat$spprich),1))
	pred <- predict(dec.mod,x.val,se.fit=T)
	seup <- pred$fit+pred$se.fit
	sedown <- pred$fit-pred$se.fit
	lines(x.val$spprich,pred$fit,col="red",lwd=2)
	lines(x.val$spprich,seup,col="red",lty=2)
	lines(x.val$spprich,sedown,col="red",lty=2)

#
# Write the plot to file
#

	dev.off()

