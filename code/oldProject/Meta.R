###	Meta-analysis of landscape declines

#
# required packages
#

	library(reshape2)	
	library(extrafont)
	library(ggplot2)
	library(maptools)
	library(classInt)

#
# Clear workspace
#

	rm(list = ls())

#
# Import and setup data
#

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

	lin <- lm(mean ~ hex, data = meta.dec)
	summary(lin)
	meta.dec$residuals <- lin$residuals

	png(filename = file.path("LRG","figs","meta","RecEffect.png")
		, width = 3000
		, height = 3000
		, res = 300
		, bg = "transparent"
		)

	par(mar = c(5, 4, 0.5, 0.5) + 0.1
		, family = "Century Gothic"
		)

	plot(meta.dec$mean ~ meta.dec$hex
		, pch = 1
		, cex = 0.6
		, xlab = "Total number of bird records (record = a species in a hexagon in a year)"
		, ylab = "Proportion of bird species declining"
		)
	
	new.x <- data.frame(
		hex = seq(min(meta.dec$hex), max(meta.dec$hex), 1)
		)

	pred.y <- predict(lin, new.x, se.fit = T)

	abline(lin, col = "red")

	lines(new.x[,1], pred.y$fit + pred.y$se.fit, lty = 2, col = "red")
	lines(new.x[,1], pred.y$fit - pred.y$se.fit, lty = 2, col = "red")

	placement <- pointLabel(x = meta.dec$hex
				, y = meta.dec$mean
				, labels = meta.dec$poly
				, cex = 0.6
				, doPlot = F
				, offset = 0
				, allowSmallOverlap = F
				)

	text(placement, labels = meta.dec$poly, cex=0.6)

	title(sub = paste("slope = "
		, round(coef(summary(lin))["hex","Estimate"],5)
		, ". p = "
		, round(coef(summary(lin))["hex","Pr(>|t|)"],5)
		, sep = "")
		)
		
	dev.off()

#
# assign action class based on level of decline
#

#
# reveg cutoffs based on mean + or - (1 * standard error)
#

	#pred.old <- predict(lin, se.fit = T)
	#meta.dec$pred <- pred.old$fit
	#meta.dec$predse <- pred.old$se.fit

	#meta.dec$cut1 <- mapply(function(x,y) x + y
	#	, x = meta.dec$pred
	#	, y = meta.dec$predse
	#	)
	
	#meta.dec$cut2 <- mapply(function(x,y) x - y
	#	, x = meta.dec$pred
	#	, y = meta.dec$predse
	#	)

	#meta.dec$decgroup <- mapply(function(x,y,z)
	#	if (is.na(z)) NA
	#	else if (z >= x) 3
	#	else if (z >= y) 2
	#	else if (z < y) 1
	#	, x = meta.dec$cut1
	#	, y = meta.dec$cut2
	#	, z = meta.dec$mean
	#	)

#
# reveg cutoffs based on natural breaks within decline residuals
#

	inter <- classIntervals(meta.dec$residuals
		, n = 4
		, style = "hclust"
		)

	meta.dec$decgroup <- sapply(X = meta.dec$residuals, FUN = function(x)
		if (is.na(x)) 0
		else if (x >= inter$brks[4]) 4
		else if (x >= inter$brks[3]) 3
		else if (x >= inter$brks[2]) 2
		else 1
		)

	meta.dec$deccolour <- sapply(X = meta.dec$decgroup, FUN = function(x)
		if (is.na(x)) 0
		else if (x == 1) "green"
		else if (x == 2) "yellow"
		else if (x == 3) "orange"
		else if (x == 4) "red"
		else 0
		)

#
# write results to file
#

	dir.create(file.path("LRG","tbl","meta"))
	write.table(meta.dec,file.path("LRG","tbl","meta","meta.csv"),row.names=F,append=F,sep=",",col.names=T)
	
#
# Print caterpillar plot
#

	dir.create(file.path("LRG","figs","meta"))

	meta.dec <- meta.dec[order(as.numeric(meta.dec$mean), decreasing = TRUE),]

	x <- meta.dec$mean
	xse <- meta.dec$se

	resolution <- 100
		
	plot.height <- (par()$mar[1] + par()$mar[2] + nrow(meta.dec)*par()$cin[2])*resolution
			
	png(file=file.path("LRG","figs","meta","meta.png")
		, width = resolution*10
		, height = plot.height
		, res = resolution
		, bg = "transparent"
		)

	par(mar = c(5, 4, 0, 1), family = "Century Gothic")
		
	dotchart(x
		, labels = meta.dec$poly
		, xlim = (c(min(x-xse), max(x+xse)))
		, xlab = "Proportion of species declining from meta analysis"
		)

	#title(main = "Summary across landscapes")
	title(sub = "Error bars are +/- standard error")

	segments(x-xse
		, 1:length(x)
		, x+xse
		, 1:length(x)
		)

	dev.off()
	graphics.off()




