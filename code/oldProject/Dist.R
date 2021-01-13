###
### Distribution
###

#
# Clear R workspace
#

	rm(list = ls())

#
# load required packages
#

	library(extrafont)		#fonts for graphics
	library(classInt)			# find 'natural' breaks in data

#
# Import and setup data
#

	bigdat <- read.csv(file.path("LRG","tbl","dataprep","raw_withfixednames.csv"))

	head(bigdat)

#
# Remove unnamed, blank (na) polygons and marine spp. Change year to before/after 2000
#

	cutoff <- 2000

	dat.remove <- which(bigdat$poly == "")
	bigdat <- bigdat[-dat.remove,]

	bigdat <- subset(bigdat, bigdat$native == 1)
	bigdat <- subset(bigdat, bigdat$type != 2)

	bigdat$poly <- factor(bigdat$poly)
	bigdat$com <- factor(bigdat$com)
	bigdat$type <- factor(bigdat$type)

	bigdat$year <- sapply(bigdat$year, function(x) if (x <cutoff) "early" else "late")
	names(bigdat)[4] <- "time"

	dat.1 <- aggregate(records ~ poly + type + com + time + hex, data = bigdat, FUN = length)

	dat.2 <- aggregate(records ~ poly + type + com + time, data = dat.1, FUN = length)
	names(dat.2)[length(names(dat.2))] <- "spphex"
	dat.2 <- subset(dat.2, spphex >2)
	
	dat.3 <- aggregate(records ~ poly + type + hex, data = dat.1, FUN = length)

	dat.4 <- aggregate(records ~ poly + type, data = dat.3, FUN = length)
	#dat.4 <- aggregate(records ~ poly + type, data = dat.4, FUN = sum)
	names(dat.4)[length(names(dat.4))] <- "lshex"

	dat.5 <- merge(dat.4,dat.2)

	dat.5$prop <- mapply(function(x,y) x/y, x = dat.5$spphex, y = dat.5$lshex)

	library(reshape2)
	dat.6 <- dcast(poly + type + com ~ time, data = dat.5)

	dat.6$dist <- mapply(function(x,y) x/y, x = dat.6$early, y = dat.6$late)
	
	dat.6 <- na.omit(dat.6)
	dat.6$poly <- factor(dat.6$poly)
	dat.6 <- data.frame(dat.6,distdec=9999)

	dat.7 <- dat.6[0,]
	
#
# split the dist data into 'natural groups' (within landscapes)
#

	for (i in 1:length(levels(dat.6$poly))){

		dat.6.sub <- dat.6[dat.6$poly == levels(dat.6$poly)[i],]
	
		inter <- classIntervals(dat.6.sub$dist
			, style = "hclust"
			)

		for (j in 1:nrow(dat.6.sub)){
			dat.6.sub$distdec[j] <- sapply(X = dat.6.sub$dist[j], FUN = function(x)
				if (x < inter$brks[2]) "black" else "grey")
				}

		dat.7 <- rbind(dat.7,dat.6.sub)
		}

#
# Ensure directories exist for ctll code outputs and empty them before running...
#

	dir.now <- file.path(file.path("LRG","tbl","dist"))
	dir.create(dir.now)
	file.list <- dir(dir.now, recursive = T)
	file.remove(file.path(dir.now,file.list))

	dir.now <- file.path("LRG","figs","dist")
	dir.create(dir.now)
	file.list <- dir(dir.now, recursive = TRUE)
	file.remove(file.path(dir.now,file.list))

	dir.create(file.path("LRG","figs","dist","LSSummary"))
	dir.create(file.path("LRG","figs","dist","SppSummary"))

#
# Export results table to .csv
#

	write.table(dat.7,file=file.path("LRG","tbl","dist","dist.csv"),row.names=F,append=F,sep=",",col.names=T)

#
# Print caterpillar plots for landscape summaries
#

	landscapes <- levels(factor(dat.7$poly))

	for (i in 1:length(landscapes)){
		dat.sub <- subset(dat.7,dat.7$poly == landscapes[i])
		dat.sub <- dat.sub[order(dat.sub$dist, decreasing = TRUE),]

		x <- dat.sub$dist

		resolution <- 100

		plot.height <- (par()$mar[1] + par()$mar[2] + nrow(dat.sub)*par()$cin[2])*resolution
			
		png(file=file.path("LRG","figs","dist","LSSummary",paste(landscapes[i],".png",sep=""))
			, width=resolution*10, height=plot.height, res=resolution
			, family = "Century Gothic"
			)
			
		label.font <- dat.sub$distdec
			
		dotchart(x
			, labels = dat.sub$com
			, xlab = "Change in distribution (pre-2000 vs. post-2000)"
			, color = label.font
			)

		title(main = landscapes[i])

		cutoff <- max(dat.sub[dat.sub[,"distdec"] == "black","dist"])
		
		abline(v=cutoff, col="red")	
	
		dev.off()
	}

#
# Print caterpillar plots for species summaries
#

	spp <- levels(factor(dat.7$com))

	for (i in 1:length(spp)){
		dat.sub <- subset(dat.7, dat.7$com == spp[i])
		dat.sub <- dat.sub[order(dat.sub$dist, decreasing = TRUE),]

		dat.sub <- na.omit(dat.sub)

		x <- dat.sub$dist

		resolution <- 100

		plot.height <- (par()$mar[1] + par()$mar[2] + nrow(dat.sub)*par()$cin[2])*resolution
			
		png(file=file.path("LRG","figs","dist","SppSummary",paste(spp[i],".png",sep=""))
			, width=resolution*10
			, height=plot.height
			, res=resolution
			, family = "Century Gothic"
			)
			
		#label.font <- sapply(dat.sub[,3], function(x) if(x < 0.05) "black" else "grey")
			
		dotchart(x
			, labels = dat.sub$poly
			#, xlim = (c(min(dat.sub[,4]-dat.sub[,5]), max(dat.sub[,4]+dat.sub[,5])))
			, xlab = "Change in distribution (pre-2000 vs. post-2000)"
			#, color = label.font
			)

		title(main = spp[i])
		#title(sub = "Landscapes with year effect p > 0.05 are grey. Error bars are +/- s.e.")
		#segments(x-dat.sub[,5],1:length(dat.sub[,4])
		#	, x+dat.sub[,5],1:length(dat.sub[,4])
		#	, col = label.font
		#	)
		
		abline(v=1, col="red")	
	
		dev.off()
	}

