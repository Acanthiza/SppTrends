###	Calculate Extinction Probability

#
# required packages
#

	library(extrafont)
	library(reshape)

#
# Clear workspace
#

rm(list = ls())

#
# Import and setup data
#

	bigdat <- read.csv(file.path("LRG","tbl","dataprep","raw_withfixednames.csv"))
	
	#rm(list=ls())
	#require(RODBC)
	#con <- odbcConnectAccess("../rel/LRG_rel.mdb")
	#dat <- sqlFetch(con,"q00140")
	#odbcCloseAll()
	#rm(con)

	head(bigdat)

	dat.remove <- which(bigdat$poly == "")
	bigdat <- bigdat[-dat.remove,]
	bigdat <- subset(bigdat, bigdat$native == 1)
	bigdat <- subset(bigdat, bigdat$type != 2)

	dat.0 <- aggregate(sites ~ poly + com, data = bigdat, FUN = sum)
	dat.0 <- subset(dat.0, dat.0$sites > 2 & dat.0$sites < 15)

	bigdat.use <- merge(bigdat,dat.0[,1:2])

	dat.1 <- aggregate(sites ~ poly + com + year + hex + type + native, data = bigdat.use, FUN = length)
	
	dat.2 <- aggregate(hex ~ poly + com + year + type + native, data = dat.1, FUN = length)

	dat.3 <- aggregate(year ~ poly + com, data = dat.2, FUN = length)
	dat.3 <- subset(dat.3, dat.3$year > 2 & dat.3$year < 15)

	dat.2 <- merge(dat.2,dat.3[,1:2])

	dat <- data.frame(poly = dat.2$poly, com = dat.2$com, year = dat.2$year, hex = dat.2$hex)

	head(dat)
	summary(dat)

	analyse.levels <- aggregate(year ~ poly + com ,dat,length)
	analyse.levels <- analyse.levels[,c(1:2)]

	analyse.levels[,1] <- factor(analyse.levels[,1])

	res <- array(NA,dim=c(length(row.names(analyse.levels)),10))
	res <- data.frame(res)
	names(res) <- c("poly","com","n","minyear","maxyear","T","Tn","McInerny","Solow","fifty")
	head(res)

	now.year <- as.numeric(format(Sys.Date(), format="%Y"))
	
	for (i in 1:nrow(res)){
		x <- as.vector(analyse.levels[i,1])
		y <- as.vector(analyse.levels[i,2])
		sub.dat <- subset(dat, poly==x & com==y)
		res[i,1] <- x
		res[i,2] <- y
		res[i,3] <- length(sub.dat$year)
		res[i,4] <- min(sub.dat$year)
		res[i,5] <- max(sub.dat$year)
		res[i,6] <- now.year-res[i,4]
		res[i,7] <- res[i,5]-res[i,4] + 1
	}

#	res <- subset(res, n>3 & Tn>2 & minyear < 2000)

	for (i in 1:nrow(res)){
		res$fifty[i] <- if (now.year - res$maxyear[i] > 50) 0 else 1
		x <- as.vector(analyse.levels[i,1])
		y <- as.vector(analyse.levels[i,2])
		sub.dat <- subset(dat, poly==x & com==y)
		res$McInerny[i] <- if (now.year - res$maxyear[i] < 2) 1.000 else (1-((res$n[i]-1)/res$Tn[i]))^(res$T[i]-res$Tn[i])
		res$Solow[i] <- if (now.year - res$maxyear[i] < 2) 1.000 else (res$Tn[i]/res$T[i])^res$n[i]
	}
		
#
# ensure directories exist for outputs from ext code, and empty them
#

	dir.now <- file.path("LRG","tbl","ext")
	dir.create(dir.now)	
	file.list <- dir(dir.now, recursive = T)
	file.remove(file.path(dir.now,file.list))

	dir.now <- file.path("LRG","figs","ext")
	dir.create(dir.now)
	file.list <- dir(dir.now, recursive = T)
	file.remove(file.path(dir.now,file.list))

	dir.create(file.path(dir.now,"LSSummary"))
	dir.create(file.path(dir.now,"SppSummary"))

#
# export results
#

	write.table(res,file=file.path("LRG","tbl","ext","ext.csv"),row.names=F,append=F,sep=",",col.names=T)
		
#
# Print caterpillar plots for landscape summaries
#

	res[,1] <- factor(res[,1])
	res[,2] <- factor(res[,2])

	for (i in 1:length(levels(res[,1]))){
		dat.sub <- subset(res, res[,1] == levels(res[,1])[i])
		dat.sub <- dat.sub[order(dat.sub[,8], decreasing = TRUE),]

		x <- dat.sub[,8]
			
		resolution <- 100

		plot.height <- (par()$mar[1] + par()$mar[2] + nrow(dat.sub)*par()$cin[2])*resolution
			
		png(file=file.path("LRG","figs","ext","LSSummary",paste(dat.sub[1,1],".png",sep=""))
			, width=resolution*10, height=plot.height, res=resolution
			, family = "Century Gothic"
			)
			
			
		label.font <- sapply(x, function(x) if(x < 0.05) "black" else "grey")
			
		dotchart(x
			, labels = dat.sub[,2]
			, xlim = c(0,1)
			, xlab = "Probability of extinction"
			, color = label.font
			)

		title(main = dat.sub[1,1])
		title(sub = "Species with ext p > 0.05 are grey.")
		
		abline(v=0.05, col="red")
		abline(v=0.1, col="orange")	

		dev.off()
	}

#
#	Print caterpillar plots for species summaries
#

	for (i in 1:length(levels(res[,2]))){
		dat.sub <- subset(res,res[,2] == levels(res[,2])[i])
		dat.sub <- dat.sub[order(dat.sub[,8], decreasing = TRUE),]

		x <- dat.sub[,8]
		
		resolution <- 100

		plot.height <- (par()$mar[1] + par()$mar[2] + nrow(dat.sub)*par()$cin[2])*resolution
			
		png(file=file.path("LRG","figs","ext","SppSummary",paste(dat.sub[1,2],".png",sep=""))
			, width=resolution*10, height=plot.height, res=resolution
			, family = "Century Gothic"
			)
				
		label.font <- sapply(x, function(x) if(x < 0.05) "black" else "grey")
			
		dotchart(x
			, labels = dat.sub[,1]
			, xlim = c(0,1)
			, xlab = "Probability of extinction"
			, color = label.font
			)

		title(main = dat.sub[1,2])
		title(sub = "Landscapes with ext p > 0.05 are grey.")
		
		abline(v=0.05, col="red")
		abline(v=0.1, col="orange")	
		
		dev.off()
	}






