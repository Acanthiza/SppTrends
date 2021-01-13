###	Decline - Calculate slope and p-value for slope


#
# Import and setup data
#

	rm(list=ls())
	require(RODBC)
	con <- odbcConnectAccess("LRG/rel/LRG_rel.mdb")
	dat <- sqlFetch(con,"q01210")
	odbcCloseAll()
	rm(con)

	head(dat)

	names(dat) <- c("poly","type","name","year","prop")

	head(dat)

#
# create directory for graph results
#

	dir.create(file.path("LRG","figs","dec"))


	library(reshape)
	poly.name <- cast(dat,poly + name + type ~ year,mean,value = "prop")
	poly.name <- poly.name[,c(1:3)]

	res <- array(NA,dim=c(length(row.names(poly.name)),9))
	res <- data.frame(res)
	names(res) <- c("run","poly","name","type","mod","n","slope","r2","p")
	head(res)

	for (i in 1:length(row.names(poly.name))){
		x <- as.vector(poly.name[i,1])
		y <- as.vector(poly.name[i,2])
		sub.dat <- subset(dat,poly==x & name==y)
		lin <- lm(prop ~ year,data = sub.dat)
		n <- nrow(sub.dat)
		m <- lin$coefficient[2]
		r2 <- summary.lm(lin)$r.squared
		p <- anova(lin)[1,5]
		res[i,] <- c(run,x,y,type=sub.dat$type[1],mod="lin",n,m,r2,p)
		minyear <- min(sub.dat$year)
		maxyear <- max(dat$year)

		
	#
	# Create graph
	#

		png(file=file.path("LRG","figs","Dec",paste(x,"_",y,".png",sep="")),width=3000,height=2000,res=200)
		

		if(p<0.1){
			plot(sub.dat$prop~sub.dat$year,xlab="Year",ylab="Proportion of Records",xlim=c(minyear,maxyear),ylim=c(0,1),sub=paste("linear model p=",round(p,3),". fitted line (red) is linear model +/- s.e.",sep=""))
			pred <- predict(lin,se.fit=T)
			seup <- pred$fit+pred$se.fit
			sedown <- pred$fit-pred$se.fit
			lines(sub.dat$year,pred$fit,col="red",lwd=2)
			lines(sub.dat$year,seup,col="red",lty=2)
			lines(sub.dat$year,sedown,col="red",lty=2)
			title(paste(x,y,sep=" - "))
		}
		
		else{
			plot(sub.dat$prop~sub.dat$year,xlab="Year",ylab="Proportion of Records",xlim=c(minyear,maxyear),ylim=c(0,1),sub=paste("linear model p=",round(p,3),sep=""))
			title(paste(x,y,sep=" - "))
		}
	
		dev.off()
	}

#
# export results
#

	write.table(res,file="LRG/tbl/dec.csv",row.names=F,append=F,sep=",",col.names=T)
