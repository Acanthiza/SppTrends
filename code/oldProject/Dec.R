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
	dat.5 <- dat.5[dat.5$spphex > 2,]										# set minimum value of 3

	dat.4 <- merge(dat.4,dat.5[,1:3])										# limit dat.4 to dat.5

	dat.7 <- merge(dat.4,dat.3)											# combine lshex and spphex

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

	dir.create(file.path("LRG","tbl","dec"))

	dir.create(file.path("LRG","figs","dec"))
	dir.now <- file.path("LRG","figs","dec")
	file.list <- dir(file.path("LRG","figs","dec"), recursive = T)
	file.remove(file.path(dir.now,file.list))

	dir.create(file.path("LRG","figs","dec","LSSummary"))
	dir.create(file.path("LRG","figs","dec","SppSummary"))

	landscapes <- factor(levels(dat$poly))
	
	for (i in 1:length(landscapes)){
		dir.create(file.path("LRG","figs","dec",landscapes[i]))
		}

#
# Create data frame as framework for subsetting and for receiving results
#

	poly.name <- cast(dat,poly + com + type ~ year,mean,value = "prop")
	poly.name <- poly.name[,c(1:3)]
	poly.name$poly <- factor(poly.name$poly)
	poly.name$com <- factor(poly.name$com)

	res <- array(NA,dim=c(length(row.names(poly.name)),9))
	res <- data.frame(res)
	names(res) <- c("poly","com","type","mod","n","slope","slopese","r2","p")
	head(res)

	for (i in 1:length(row.names(poly.name))){
		x <- as.vector(poly.name[i,1])
		y <- as.vector(poly.name[i,2])
		sub.dat <- subset(dat,poly==x & com==y)
		lin <- lm(prop ~ year,data = sub.dat)
		n <- nrow(sub.dat)
		m <- lin$coefficient[2]
		m.se <- coef(summary(lin))["year","Std. Error"]
		r2 <- summary.lm(lin)$r.squared
		p <- anova(lin)[1,5]
		res[i,] <- c(x,y,type=sub.dat$type[1],mod="lin",n,m,m.se,r2,p)
		minyear <- min(sub.dat$year)
		maxyear <- max(dat$year)
		
	#
	# Create graph
	#

		png(file=file.path("LRG","figs","dec",x,paste(y,".png",sep=""))
			,width=3000,height=2000,res=200
			,family = "Century Gothic"
			)
		

		if(p<0.1){
			plot(sub.dat$prop~sub.dat$year
				, xlab="Year"
				, ylab="Proportion of Records"
				, xlim=c(minyear,maxyear)
				, ylim=c(0,1)
				, sub=paste("linear model p=",round(p,3),". fitted line (red) is linear model +/- s.e.",sep="")
				)
			x.val <- data.frame(year=seq(minyear,maxyear,0.1))
			pred <- predict.lm(lin,x.val,type="resp",se=T)
			lines(x.val[,1],pred$fit,col="red",lwd=2)
			lines(x.val[,1],pred$fit+pred$se.fit,col="red",lty=2)
			lines(x.val[,1],pred$fit-pred$se.fit,col="red",lty=2)
			title(paste(x,y,sep=" - "))
		} else {
			plot(sub.dat$prop~sub.dat$year
				, xlab="Year"
				, ylab="Proportion of Records"
				, xlim=c(minyear,maxyear)
				, ylim=c(0,1)
				, sub=paste("linear model p=",round(p,3),sep="")
				)
			title(paste(x,y,sep=" - "))
		}
	
		dev.off()
	}

#
# export results
#

	write.table(res,file.path("LRG","tbl","dec","dec.csv"),row.names=F,append=F,sep=",",col.names=T)

#
# Print caterpillar plots for landscape summaries
#

	landscapes <- factor(levels(poly.name[,1]))

	for (i in 1:length(landscapes)){
		dat.sub <- subset(res, res$poly == landscapes[i])
		dat.sub <- dat.sub[order(as.numeric(dat.sub$slope), decreasing = TRUE),]

		x <- as.numeric(dat.sub$slope)
		xse <- as.numeric(dat.sub$slopese)

		resolution <- 100

		plot.height <- (par()$mar[1] + par()$mar[2] + nrow(dat.sub)*par()$cin[2])*resolution
			
		png(file=file.path("LRG","figs","dec","LSSummary",paste(landscapes[i],".png",sep=""))
			, width=resolution*10, height=plot.height, res=resolution
			, family = "Century Gothic"
			)
			
		label.font <- sapply(dat.sub$p, function(x) if(x < 0.05) "black" else "grey")
			
		dotchart(x
			, labels = dat.sub$com
			, xlim = (c(min(x-xse), max(x+xse)))
			, xlab = "Effect of year on presence"
			, color = label.font
			)

		title(main = dat.sub[1,1])
		title(sub = "Species with year effect p > 0.05 are grey. Error bars are +/- s.e.")
		segments(x-xse,1:length(dat.sub[,4])
			, x+xse,1:length(dat.sub[,4])
			, col = label.font
			)

		abline(v=0, col="red")	
	
		dev.off()
	}

#
# Print caterpillar plots for species summaries
#

	spp <- factor(levels(poly.name[,2]))

	for (i in 1:length(spp)){
		dat.sub <- subset(res, res$com == spp[i])
		dat.sub <- dat.sub[order(as.numeric(dat.sub$slope), decreasing = TRUE),]

		x <- as.numeric(dat.sub$slope)
		xse <- as.numeric(dat.sub$slopese)
		
		resolution <- 100

		plot.height <- (par()$mar[1] + par()$mar[2] + nrow(dat.sub)*par()$cin[2])*resolution
			
		png(file=file.path("LRG","figs","dec","SppSummary",paste(spp[i],".png",sep=""))
			, width=resolution*10, height=plot.height, res=resolution
			, family = "Century Gothic"
			)

		label.font <- sapply(dat.sub$p, function(x) if(x < 0.05) "black" else "grey")
			
		dotchart(x
			, labels = dat.sub[,1]
			, xlim = (c(min(x-xse), max(x+xse)))
			, xlab = "Effect of year on presence"
			, color = label.font
			)

		title(main = spp[i])
		title(sub = "Landscapes with year effect p > 0.05 are grey. Error bars are +/- s.e.")
		segments(x-xse,1:length(dat.sub[,4])
			, x+xse,1:length(dat.sub[,4])
			, col = label.font
			)
		abline(v=0, col="red")	
	
		dev.off()
	}


