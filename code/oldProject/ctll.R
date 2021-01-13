###
### List length
###

#
# clear workspace
#

	rm(list=ls())


#
# Increase memory size
#

	#memory.limit(2048)

#
# load required libraries
#

	library(extrafont)

#
# Import and setup data
#

	bigdat <- read.csv(file.path("LRG","tbl","dataprep","raw_withfixednames.csv"))

	#rm(list=ls())
	#require(RODBC)
	#con <- odbcConnectAccess("LRG/Rel/LRG_rel.mdb")
	#dat <- sqlFetch(con,"q01335")
	#odbcCloseAll()
	#rm(con)

	head(bigdat)

#
# Remove unnamed, blank (na) polygons
#

	dat.remove <- which(bigdat$poly == "")
	bigdat <- bigdat[-dat.remove,]
	bigdat <- subset(bigdat, native == 1)
	bigdat <- subset(bigdat, type != 2)
	bigdat$list <- paste(bigdat$hex,bigdat$year,bigdat$type)

	dat.1 <- aggregate(records ~ poly + type + com + year + hex + list, data = bigdat, FUN = length)

	dat <- data.frame(poly = dat.1$poly
		, hex = dat.1$hex
		, year = dat.1$year
		, list = dat.1$list
		, type = dat.1$type
		, com = dat.1$com
		, n = dat.1$records
		)
	
	dat.check <- aggregate(n ~ year, data = dat, FUN = length)
	plot(dat.check)	# clear increase in records from 1982 onwards

#
# Set earliest year to 1982 based on previous plot
#

	dat <- subset(dat, dat$year > 1981)

	listlength <- aggregate(n ~ poly + year + list + type, data = dat, FUN = length)
	head(listlength)
	names(listlength) <- c("poly","year","list","type","ll")
	listlength <- subset(listlength, ll > 1)

#
# Set minimum number of lists per landscape a species appears on
#

	ls.spplist <- aggregate(n ~ poly + type + com, data = dat, FUN = length)
	head(ls.spplist)
	names(ls.spplist) <- c("poly","type","com","records")
	ls.spplist <- subset(ls.spplist, ls.spplist$records > 9)

	ls.poss.lists <- merge(ls.spplist[,1:3],listlength)
	head(ls.poss.lists)

	ls.lists <- merge(ls.poss.lists,dat,all.x = T)
	ls.lists <- ls.lists [order(ls.lists$list, decreasing = TRUE),]
	ls.lists[1:200,]

	dat <- data.frame(poly=ls.lists$poly
		,time=ls.lists$year
		,type=ls.lists$type
		,list=ls.lists$list
		,ll=ls.lists$ll
		,com=ls.lists$com
		,p=ls.lists$n
		)

#
# Replace na with 0 in presence field and set up absence field
#

	dat[,7][is.na(dat[,7])]<-0
		
	dat$a <- sapply(dat$p, function(x) if(x == 0) 1 else 0)

#
# Minimum number of years per species (within landscape and type)
#

	dat.1 <- aggregate(p ~ poly + time + type + com, data = dat, FUN = length)
	dat.2 <- aggregate(time ~ poly + type + com, data = dat.1, FUN = length)
	dat.3 <- subset(dat.2, time > 4)
	head(dat.3)

	dat <- merge(dat.3[,1:3],dat)

#
# Create data frame with levels to run through analyses
#

	analyse.levels <- aggregate(p ~ poly + com, data = dat, length)
	analyse.levels$year.coeff <- 1
	analyse.levels$year.se <- 1
	
	analyse.levels <- analyse.levels [ order(analyse.levels [,1], analyse.levels [,2]),]

	analyse.levels[,1] <- factor(analyse.levels[,1])
	analyse.levels[,2] <- factor(analyse.levels[,2])

	head(analyse.levels)

#
# Ensure directories exist for ctll code outputs and empty them before running...
#

	dir.create(file.path("LRG","tbl","ctll"))

	dir.create(file.path("LRG","figs","ctll"))
	dir.now <- (file.path("LRG","figs","ctll"))
	file.list <- dir(file.path("LRG","figs","ctll"), recursive = TRUE)
	file.remove(file.path(dir.now,file.list))

	dir.create(file.path("LRG","figs","ctll","LSSummary"))
	dir.create(file.path("LRG","figs","ctll","SppSummary"))
	for (i in 1:length(levels(analyse.levels[,1]))){
		dir.create(file.path("LRG","figs","ctll",levels(analyse.levels[,1])[i]))
		}

#
# Run model
#
		
	for (i in 1:nrow(analyse.levels)){
		x <- analyse.levels[i,1]
		y <- analyse.levels[i,2]
		dat.sub <- subset(dat, as.character(dat$poly) == x)
		dat.sub <- subset(dat.sub, as.character(dat.sub$com) == y)

		bm <- glm(cbind(p, a) ~ ll + time, family=binomial(logit), data=dat.sub, maxit = 50)

		analyse.levels[i,3] <- coef(summary(bm))["time","Pr(>|z|)"]
		analyse.levels[i,4] <- coef(summary(bm))["time","Estimate"]
		analyse.levels[i,5] <- coef(summary(bm))["time","Std. Error"]
		x.val <- data.frame(time=seq(min(dat.sub$time),max(dat.sub$time),0.1))
		ct.ll <- mean(dat.sub$ll)
		x.val <- data.frame(ll=ct.ll,x.val)
		pred.bm <- predict.glm(bm,x.val,type="resp",se=T)
		
		png(file=file.path("LRG","figs","ctll",x,paste(y,".png",sep=""))
			, width=2000, height=2000, res=200
			, family = "Century Gothic"
			, bg = "transparent"
			)	
		
		plot(dat.sub$time
			, dat.sub$p
			, pch=""
			#, cex=0.5
			, col="red"
			, ylab="P[obs]"
			, xlab = "Year"
			, xlim=c(min(dat.sub$time)
			, max(dat.sub$time))
			, ylim=c(0,1)
			)
			
		lines(x.val$time, pred.bm$fit, type="l", col="red")
		lines(x.val$time, pred.bm$fit - pred.bm$se.fit, lty=2, col="red")
		lines(x.val$time, pred.bm$fit + pred.bm$se.fit, lty=2, col="red")		

		title(main = paste(x,y, sep = ": "))
		title(sub = paste("At mean list length = ",round(ct.ll,1),sep=""))

		dev.off()	
	}

	write.table(analyse.levels,file=file.path("LRG","tbl","ctll","ctll.csv"),row.names=F,append=F,sep=",",col.names=T)

	head(analyse.levels)

	#analyse.levels <- read.csv(file=file.path("LRG","tbl","ctll","ctll.csv"))

#
# Print caterpillar plots for landscape summaries
#

	for (i in 1:length(levels(analyse.levels[,1]))){
		dat.sub <- subset(analyse.levels,analyse.levels[,1] == levels(analyse.levels[,1])[i])
		dat.sub <- dat.sub[order(dat.sub[,4], decreasing = TRUE),]

		#xval <- mapply(function(x,y) ifelse(x < 0.05, y, NA), x = dat.sub[,3], y = dat.sub[,4])

		x <- dat.sub[,4]
		#x <- xval			

		resolution <- 100

		plot.height <- (par()$mar[1] + par()$mar[2] + nrow(dat.sub)*par()$cin[2])*resolution
			
		png(file=file.path("LRG","figs","ctll","LSSummary",paste(dat.sub[1,1],".png",sep=""))
			, width=resolution*10, height=plot.height, res=resolution
			, family = "Century Gothic"
			)
			
		label.font <- sapply(dat.sub[,3], function(x) if(x < 0.05) "black" else "grey")
			
		dotchart(x
			, labels = dat.sub[,2]
			#, xlim = (c(min(dat.sub[,4]-dat.sub[,5]), max(dat.sub[,4]+dat.sub[,5])))
			, xlim = c(-1.5,1.5)
			, xlab = "Effect of year on presence"
			, color = label.font
			)

		title(main = dat.sub[1,1])
		title(sub = "Species with year effect p > 0.05 are grey. Error bars are +/- s.e.")
		segments(x-dat.sub[,5],1:length(dat.sub[,4])
			, x+dat.sub[,5],1:length(dat.sub[,4])
			, col = label.font
			)

		abline(v=0, col="red")	
	
		dev.off()
	}

#
# Print caterpillar plots for species summaries
#

	for (i in 1:length(levels(analyse.levels[,2]))){
		dat.sub <- subset(analyse.levels,analyse.levels[,2] == levels(analyse.levels[,2])[i])
		dat.sub <- dat.sub[order(dat.sub[,4], decreasing = TRUE),]

		#xval <- mapply(function(x,y) ifelse(x < 0.05, y, NA), x = dat.sub[,3], y = dat.sub[,4])

		x <- dat.sub[,4]
		#x <- xval			

		resolution <- 100

		plot.height <- (par()$mar[1] + par()$mar[2] + nrow(dat.sub)*par()$cin[2])*resolution
			
		png(file=file.path("LRG","figs","ctll","SppSummary",paste(dat.sub[1,2],".png",sep=""))
			, width=resolution*10, height=plot.height, res=resolution
			, family = "Century Gothic"
			)
			
			
		label.font <- sapply(dat.sub[,3], function(x) if(x < 0.05) "black" else "grey")
			
		dotchart(x
			, labels = dat.sub[,1]
			#, xlim = (c(min(dat.sub[,4]-dat.sub[,5]), max(dat.sub[,4]+dat.sub[,5])))
			, xlim = c(-1.5,1.5)
			, xlab = "Effect of year on presence"
			, color = label.font
			)

		title(main = dat.sub[1,2])
		title(sub = "Landscapes with year effect p > 0.05 are grey. Error bars are +/- s.e.")
		segments(x-dat.sub[,5],1:length(dat.sub[,4])
			, x+dat.sub[,5],1:length(dat.sub[,4])
			, col = label.font
			)
		abline(v=0, col="red")	
	
		dev.off()
	}


