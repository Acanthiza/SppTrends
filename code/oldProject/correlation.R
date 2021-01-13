###	Decline - Calculate slope and p-value for slope

#
# Clear workspace
#

	rm(list=ls())

#
# load required libraries
#

	library(extrafont)
	library(reshape)
	library(Hmisc)

#
# Import and setup data
#

	bigdat <- read.csv(file.path("LRG","tbl","dataprep","raw_withfixednames.csv"))

	head(bigdat)
		
	bigdat$hex <- as.factor(bigdat$hex)
	bigdat$type <- as.factor(bigdat$type)
	bigdat$native <- as.factor(bigdat$native)
	
	summary(bigdat)

	dat.remove <- which(bigdat$poly == "")
	bigdat <- bigdat[-dat.remove,]
	bigdat <- subset(bigdat, native == 1)
	bigdat <- subset(bigdat, type != 2)
	bigdat <- bigdat[bigdat$year > 1979,]
	bigdat$poly <- factor(bigdat$poly)

	dat.1 <- aggregate(records ~ poly + hex + com + year + type, data = bigdat, FUN = length)	# basic summary
	
	dat.2 <- aggregate(records ~ poly + hex + year + type, data = dat.1, FUN = length)		# n spp per hex per year
	
	dat.3 <- aggregate(records ~ poly + year + type, data = dat.2, FUN = length)			# n hex per year
	names(dat.3) <- c("poly","year","type","lshex")
	dat.3 <- dat.3[dat.3$lshex > 2, ]

	dat.1 <- merge(dat.1,dat.3[,1:3])										# limit dat.1 to dat.3

	dat.4 <- aggregate(records ~ poly + com + year + type, data = dat.1, FUN = length)		# n hex per year per spp
	names(dat.4) <- c("poly","com","year","type","spphex")

	dat.5 <- aggregate(spphex ~ poly + com + type, data = dat.4, FUN = length)			# n year per spp
	dat.5 <- dat.5[dat.5$spphex > 4,]										# set minimum value of 5 (needed for spearman) -> LS spp list

	dat.4 <- merge(dat.4,dat.5[,1:3])										# limit dat.4 to dat.5

	#dat.6 <- merge(dat.5[,1:3],dat.3)										# combine lshex and LS spp list

	dat.7 <- merge(dat.3,dat.4,all.x = T)									# add spphex to dat.6 or dat.3

	dat.7[is.na(dat.7)] <- 0											# replace NA from last merge with 0 (no records)

	dat.7$prop <- mapply(function(x,y) x/y, x = dat.7$spphex, y = dat.7$lshex)
	
	dat <- data.frame(poly = dat.7$poly
		, type = dat.7$type
		, com = dat.7$com
		, year = dat.7$year
		, prop = dat.7$prop
		)

#
# ensure directories exist for outputs from dec code
#

	dir.create(file.path("LRG","tbl","cor"))

	dir.create(file.path("LRG","figs","cor"))
	dir.now <- file.path("LRG","figs","cor")
	file.list <- dir(file.path("LRG","figs","cor"), recursive = T)
	file.remove(file.path(dir.now,file.list))

	dir.create(file.path("LRG","figs","cor","LSSummary"))
	dir.create(file.path("LRG","figs","cor","SppSummary"))

	landscapes <- factor(levels(dat$poly))
	
	for (i in 1:length(landscapes)){
		dir.create(file.path("LRG","figs","cor",landscapes[i]))
		}

#
# Create data frame as framework for subsetting and for receiving results
#

	poly.name <- dat.5[,1:3]

	res <- array(NA,dim=c(length(row.names(poly.name)),7))
	res <- data.frame(res)
	names(res) <- c("poly","com","type","mod","n","corr","p")
	head(res)

	cor.type <- "pearson"
	
	for (i in 1:length(row.names(poly.name))){
		x <- as.vector(poly.name[i,1])
		y <- as.vector(poly.name[i,2])
		sub.dat <- subset(dat,poly==x & com==y)

		corr <- rcorr(cbind(sub.dat$year,sub.dat$prop), type=cor.type) 

		n <- corr$n[1,1]

		m <- corr$r[1,2]
		
		p <- corr$P[1,2]

		res[i,] <- c(x,y,type=sub.dat$type[1],mod=cor.type,n,m,p)
		minyear <- min(sub.dat$year)
		maxyear <- max(sub.dat$year)
		
	#
	# Create graph
	#

		png(file=file.path("LRG","figs","cor",x,paste(y,".png",sep=""))
			, width=2000,height=2000,res=200
			, family = "Century Gothic"
			, bg = "transparent"
			)
		
		plot(sub.dat$prop~sub.dat$year
			, xlab="Year"
			, ylab="Recording Frequency"
			, xlim=c(minyear,maxyear)
			, ylim=c(0,1)
			)
	
		abline(lm(sub.dat$prop~sub.dat$year))

		title(main = paste(x,y,sep=": ")
			, sub = bquote(paste("Pearson's "
				, italic(r)
				, " = "
				, .(round(m, 2))
				, " with P = "
				, .(round(p, 2))
				))
			)
			
		dev.off()
	}

#
# export results
#

	write.table(res,file.path("LRG","tbl","cor","cor.csv"),row.names=F,append=F,sep=",",col.names=T)

#
# Print caterpillar plots for landscape summaries
#

	landscapes <- factor(levels(poly.name[,1]))

	for (i in 1:length(landscapes)){
		dat.sub <- subset(res, res$poly == landscapes[i])
		dat.sub <- dat.sub[order(as.numeric(dat.sub$corr), decreasing = TRUE),]

		x <- as.numeric(dat.sub$corr)

		resolution <- 100

		plot.height <- (par()$mar[1] + par()$mar[2] + nrow(dat.sub)*par()$cin[2])*resolution
			
		png(file=file.path("LRG","figs","cor","LSSummary",paste(landscapes[i],".png",sep=""))
			, width=resolution*10, height=plot.height, res=resolution
			, family = "Century Gothic"
			)
			
		label.font <- sapply(dat.sub$p, function(x) if(x < 0.05) "black" else "grey")
			
		dotchart(x
			, labels = dat.sub$com
			, xlab = "Correlation of frequency and year"
			, color = label.font
			)

		title(main = dat.sub[1,1])
		title(sub = "Species with correlation P > 0.05 are grey.")
		
		abline(v=0, col="red")	
	
		dev.off()
	}

#
# Print caterpillar plots for species summaries
#

	spp <- levels(factor(poly.name[,2]))

	for (i in 1:length(spp)){
		dat.sub <- res[res$com == spp[i],]
		dat.sub <- dat.sub[order(as.numeric(dat.sub$corr), decreasing = TRUE),]

		x <- as.numeric(dat.sub$corr)
				
		resolution <- 100

		plot.height <- (par()$mar[1] + par()$mar[2] + nrow(dat.sub)*par()$cin[2])*resolution
			
		png(file=file.path("LRG","figs","cor","SppSummary",paste(spp[i],".png",sep=""))
			, width=resolution*10, height=plot.height, res=resolution
			, family = "Century Gothic"
			, bg = "transparent"
			)

		label.font <- sapply(dat.sub$p, function(x) if(x < 0.05) "black" else "grey")
			
		dotchart(x
			, labels = dat.sub[,1]
			, xlab = "Correlation of frequency and year"
			, color = label.font
			)

		title(main = spp[i])
		title(sub = "Landscapes with correlation P > 0.05 are grey")
		
		abline(v=0, col="red")	
	
		dev.off()
	}


