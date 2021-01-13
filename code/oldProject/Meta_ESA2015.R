###	Meta-analysis of landscape declines

#
# required packages
#

	library(reshape2)	
	library(extrafont)
	library(ggplot2)
	library(maptools)
	library(classInt)
	library(plotrix)

#
# Clear workspace
#

	rm(list = ls())
	
#
# Set working directory
#

	setwd("//env/IST/SRC/EcoAnalysis/Projects/CC_Stream2_CCLA/Workspace_NW_Cbackup")
	
#
# Import and setup data
#

	ls <- read.csv("C:/Workspace_NW/ESA2015/tbl/Landscape.csv")
	
	spppoly <- read.csv(file.path("LRG","tbl","dataprep","spppoly.csv"))
	spppoly <- spppoly[spppoly$type == 1 | spppoly$type == 3,]
	
	poly.rec <- aggregate(hex ~ poly, data = spppoly, FUN = length)

	levels(spppoly$poly)[1] <- NA
	poly <- data.frame(poly=levels(spppoly$poly))
	polySR <- aggregate(com ~ poly + type, data = spppoly, FUN = length)
		
	ctll <- read.csv(file.path("LRG","tbl","ctll","ctll.csv"))
	ctll <- merge(ctll,spppoly[,1:3])
	dec <- read.csv(file.path("LRG","tbl","dec","dec.csv"))
	dist <- read.csv(file.path("LRG","tbl","dist","dist.csv"))
	ext <- read.csv(file.path("LRG","tbl","ext","ext.csv"))
	ext <- merge(ext,spppoly[,1:3])
	corr <- read.csv(file.path("LRG","tbl","cor","cor.csv"))

	ctll.sub <- ctll[ctll$p < 0.05 & ctll$year.coeff < 0,]
	ctll.dec <- aggregate(p ~ poly + type, data = ctll.sub, FUN = length)

	dec.sub <- dec[dec$p < 0.05 & dec$slope < 0,]
	dec.dec <- aggregate(p ~ poly + type, data = dec.sub, FUN = length)

	dist.sub <- dist[dist$distdec == "black",]
	dist.dec <- aggregate(dist ~ poly + type, data = dist.sub, FUN = length)

	ext.sub <- ext[ext$McInerny <0.05 | ext$fifty == 0,]
	ext.dec <- aggregate(McInerny ~ poly + type, data = ext.sub, FUN = length)

	corr.sub <- corr[corr$p < 0.05 & corr$corr < 0,]
	corr.dec <- aggregate(corr ~ poly + type, data = corr.sub, FUN = length)

	meta <- merge(polySR, ctll.dec, all.x = T)
	names(meta)[length(names(meta))] <- "ctll"
	
	meta <- merge(meta, dec.dec, all.x = T)
	names(meta)[length(names(meta))] <- "dec"

	meta <- merge(meta, dist.dec, all.x = T)
	names(meta)[length(names(meta))] <- "dist"

	meta <- merge(meta, ext.dec, all.x = T)
	names(meta)[length(names(meta))] <- "ext"

	meta <- merge(meta, corr.dec, all.x = T)
	names(meta)[length(names(meta))] <- "corr"

	meta.prop <- meta
	
	meta.prop[,4] <- mapply (function(x,y) x/y, x = meta[,4], y = meta[,3])
	meta.prop[,5] <- mapply (function(x,y) x/y, x = meta[,5], y = meta[,3])		
	meta.prop[,6] <- mapply (function(x,y) x/y, x = meta[,6], y = meta[,3])		
	meta.prop[,7] <- mapply (function(x,y) x/y, x = meta[,7], y = meta[,3])

	meta.melt <- melt(meta, id = names(meta)[1:3], na.rm = T)

	meta.melt$prop <- mapply (function(x,y) x/y, x = meta.melt$value, y = meta.melt$com)

	#meta.melt <- meta.melt[meta.melt[,2] == 1
	#	| meta.melt[,2] == 3
	#	,]

	meta.melt <- meta.melt[meta.melt$variable != "ext",]
	meta.dec.mean <- aggregate(prop ~ poly, data = meta.melt, FUN = mean)
	meta.dec.length <- aggregate(prop ~ poly, data = meta.melt, FUN = length)
	meta.dec.var <- aggregate(prop ~ poly, data = meta.melt, FUN = var)
	meta.dec <- data.frame(meta.dec.mean, meta.dec.length[,2], meta.dec.var[,2])
	names(meta.dec) <- c("poly","mean","n","var")
	meta.dec$se <- mapply(function(x,y) sqrt(x/y), x = meta.dec$var, y = meta.dec$n)

#
# test effect of sample size on decline
#

	meta.dec <- merge(meta.dec,poly.rec)
	
	meta.dec <- merge(meta.dec,ls[,1:2])

	lin <- lm(mean ~ hex, data = meta.dec)
	summary(lin)
	meta.dec$residuals <- lin$residuals
	
	meta.res <- resid(lin)
	
	
	png(filename = file.path("LRG","figs","meta","RecEffect_ESA2015.png")
		, width = 3000
		, height = 2000
		, res = 320
		, bg = "transparent"
		)

	par(mar = c(5, 4, 2, 0.5) + 0.1
		, family = "Century Gothic"
		)

	plot(meta.res ~ meta.dec$hex
		, pch = 16
		, cex.lab = 1.4
		#, xlab = "Total Bird Records"
		#, ylab = "Adaptive Capacity based on Species Decline"
		, xlab = ""
		, ylab = ""
		, yaxt = "n"
		, xaxt = "n"
		, col = color.scale(meta.res,c(0,1,1),c(1,1,0),0)
		)
	
	abline(h=0)
	
# 	mtext("residuals from regression of total records against proportion of bird species declining"
# 	      , side = 2
# 	      , cex = 0.9
# 	      , outer = T
# 	      )
	
	segments(meta.dec$hex,0,meta.dec$hex,meta.res
	         , col = color.scale(meta.res,c(0,1,1),c(1,1,0),0)
	         , lwd = 2
	)
	
# 	mtext(c("less","more"),side = 1
# 	      , at = c(min(meta.dec$hex)+((max(meta.dec$hex)-min(meta.dec$hex))*0.05)
# 	               , max(meta.dec$hex)-((max(meta.dec$hex)-min(meta.dec$hex))*0.05))
# 	      , cex = 1.2
# 	      )
	
# 	mtext(c("less","more")
# 	      , side = 2
# 	      , at = c(min(meta.res)+((max(meta.res)-min(meta.res))*0.05)
# 	               , max(meta.res)-((max(meta.res)-min(meta.res))*0.05))
# 	      , cex = 1.2
# 	      , las = 2
# 	)
	
	placement <- pointLabel(x = meta.dec$hex
				, y = meta.res
				, labels = meta.dec$ShortName
				#, cex = 0.6
				, doPlot = F
				, offset = 0
				, allowSmallOverlap = F
				)

	text(placement, labels = meta.dec$ShortName
	     #, cex=0.6
	     )

	title(main = "Residuals from regression of proportion species declining against total bird records"
	      , sub = bquote(F[list(1,33)]*.(paste(" = "
	                               , round(summary(lin)$fstatistic[1],2)
	                               , "; p = "
	                               , round(coef(summary(lin))["hex","Pr(>|t|)"],5)
	                               , "."
	                               , sep = ""))
	                     ~R^2*.(paste(" = "
	                               ,round(summary(lin)$r.squared, 3)
	                               , sep=""))
	                     )
	      , cex.sub = 0.8
	      )
	
	dev.off()

