

###	Decline - Calculate slope and p-value for slope

##Import and setup data
rm(list=ls())
require(RODBC)
con <- odbcConnectAccess("C:/Workspace_NW/LRG/rel/LRG_rel.mdb")
dat <- sqlFetch(con,"q01210")
maxRun <- sqlFetch(con,"q00090")
odbcCloseAll()
rm(con)

head(dat)

names(dat) <- c("poly","type","name","year","prop")

head(dat)

run <- maxRun + 1

## remove any existing outputs
#do.call(file.remove,list(list.files("C:/Workspace_NW/LRG/R/dec",full.names=TRUE)))

##	create directory for graph results
dir.create(paste("C:/Workspace_NW/LRG/R/dec/run",run,sep=""))

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
	if(p<0.1)
		{plot(sub.dat$prop~sub.dat$year,xlab="Year",ylab="Proportion of Records",xlim=c(minyear,maxyear),ylim=c(0,1),sub=paste("linear model p=",round(p,3),". fitted line (red) is linear model +/- s.e.",sep=""))
		pred <- predict(lin,se.fit=T)
		seup <- pred$fit+pred$se.fit
		sedown <- pred$fit-pred$se.fit
		lines(sub.dat$year,pred$fit,col="red",lwd=2)
		lines(sub.dat$year,seup,col="red",lty=2)
		lines(sub.dat$year,sedown,col="red",lty=2)
		title(paste(x,y,sep=" - "))
	}
	else
		{plot(sub.dat$prop~sub.dat$year,xlab="Year",ylab="Proportion of Records",xlim=c(minyear,maxyear),ylim=c(0,1),sub=paste("linear model p=",round(p,3),sep=""))
		title(paste(x,y,sep=" - "))
	}
	filename <- paste("C:/Workspace_NW/LRG/R/dec/run",run,"/",x,"_",y,sep="")
	savePlot(filename,"jpg")
}


# export results
write.table(res,file="C:/Workspace_NW/LRG/R/dec.csv",row.names=F,append=T,sep=",",col.names=F)


### Paste out the distribution query...
## Import and setup data
rm(list=ls())
require(RODBC)
con <- odbcConnectAccess("C:/Workspace_NW/LRG/rel/LRG_rel.mdb")
dat <- sqlFetch(con,"q01320")
maxRun <- sqlFetch(con,"q00090")
odbcCloseAll()
rm(con)

run <- maxRun + 1
dat <- cbind(run,dat)
names(dat) <- c("run","ls","com","propA","propB","propC")

# export results
write.table(dat,file="C:/Workspace_NW/LRG/R/dist.csv",row.names=F,append=T,sep=",",col.names=F)


### 	List length
setwd("C:/Workspace_NW/LRG/R")

##	Increase memory size
memory.limit(2048)

##	Import and setup data
rm(list=ls())
require(RODBC)
con <- odbcConnectAccess("C:/Workspace_NW/LRG/Rel/LRG_rel.mdb")
dat <- sqlFetch(con,"q01390")
maxRun <- sqlFetch(con,"q00090")
odbcCloseAll()
rm(con)

head(dat)

names(dat) <- c("poly","time","type","list","ll","name","p")

dat$p[is.na(dat$p)]<-0

dat$time <- as.factor(dat$time)

dat$a <- 1-dat$p

run <- maxRun + 1

## 	Remove any existing outputs
#do.call(file.remove,list(list.files("C:/Workspace_NW/LRG/R/ll",full.names=TRUE)))

##	create directory for graph results
dir.create(paste("C:/Workspace_NW/LRG/R/ll/run",run,sep=""))

##	prepare data and results data.frame
library(reshape)
poly.name <- cast(dat,poly + type + name ~ time,mean,value = "p")
poly.name <- poly.name[,c(1:3)]

res.1 <- array(NA,dim=c(length(row.names(poly.name)),15))
res.1 <- data.frame(res.1)
names(res.1) <- c("run","poly","name","type","time","ct","ctll","ctllPOcc","se","intercept","intP","slope","slopeP","PInc","PIncSE")
head(res.1)

res.2 <- array(NA,dim=c(length(row.names(poly.name)),15))
res.2 <- data.frame(res.2)
names(res.2) <- c("run","poly","name","type","time","ct","ctll","ctllPOcc","se","intercept","intP","slope","slopeP","PInc","PIncSE")

#	Define and calculate a measure of central tendancy (ct) to use in predicted values of POcc
ct <- "mean"
ls.type.ct.ll <- aggregate(dat$ll, list(poly=dat$poly,type=dat$type), FUN = mean)
names(ls.type.ct.ll) <- c("poly","type","ctll")

head(res.2)

for (i in 1:length(row.names(poly.name))){
	x <- as.vector(poly.name[i,1])
	y <- as.vector(poly.name[i,3])
	sub.dat <- subset(dat,poly==x & name==y)
	z <- as.factor(poly.name[i,2])
	
	dat.1 <- subset(sub.dat,sub.dat$time==levels(sub.dat$time)[1])
	dat.2 <- subset(sub.dat,sub.dat$time==levels(sub.dat$time)[2])
	
	bm.1 <- glm(cbind(p, a) ~ ll, family=binomial(logit), data=dat.1, maxit = 50)
	x.val <- data.frame(ll=seq(min(dat$ll),max(dat$ll),1))
	pred.bm.1 <- predict.glm(bm.1,x.val,type="resp",se=T)
	
	bm.2 <- glm(cbind(p, a) ~ ll, family=binomial(logit), data=dat.2, maxit = 50)
	pred.bm.2 <- predict.glm(bm.2,x.val,type="resp",se=T)

	plot(dat.1$ll, dat.1$p, pch=3, cex=0.5, col="red", xlab="List Length",ylab="P[obs]",xlim=c(min(sub.dat$ll),max(sub.dat$ll)),ylim=c(0,1))
	title(paste(x,y,sep=" - "))
	title(sub="Red is for lists < year 2000 and Green is >=2000. Lines are fitted logisitc regression +/- s.e.",cex=0.5)
	
	lines(x.val$ll, pred.bm.1$fit, type="l", col="red")
	lines(x.val$ll, pred.bm.1$fit - pred.bm.1$se.fit, lty=2, col="red")
	lines(x.val$ll, pred.bm.1$fit + pred.bm.1$se.fit, lty=2, col="red")
	
	points(dat.2$ll, dat.2$p, pch=4, cex=0.5, col="green")
	lines(x.val$ll, pred.bm.2$fit, type="l", col="green")
	lines(x.val$ll, pred.bm.2$fit - pred.bm.2$se.fit, lty=2, col="green")
	lines(x.val$ll, pred.bm.2$fit + pred.bm.2$se.fit, lty=2, col="green")

	ct.ll <- data.frame(ll=ls.type.ct.ll[ls.type.ct.ll$poly==x & ls.type.ct.ll$type==z,3])
	predict.ct.ll.1 <- predict.glm(bm.1,ct.ll,type="resp", se=T)
	predict.ct.ll.2 <- predict.glm(bm.2,ct.ll,type="resp", se=T)

	abline(v=ct.ll,lty=3)
	mtext(paste(ct,"list length",sep=" "), side=3, at = ct.ll)

	filename <- paste("C:/Workspace_NW/LRG/R/ll/run",run,"/",x,"_",y,sep="")
	savePlot(filename,"jpg")

	#	p-value intercept
	bm.1.int.p <- coef(summary(bm.1))["(Intercept)","Pr(>|z|)"]
	bm.2.int.p <- coef(summary(bm.2))["(Intercept)","Pr(>|z|)"]

	#	p-value slope
	bm.1.slope.p <- coef(summary(bm.1))["ll","Pr(>|z|)"]
	bm.2.slope.p <- coef(summary(bm.2))["ll","Pr(>|z|)"]

	#	increased odds - estimate (slope)
	PInc.1 <- exp(coef(summary(bm.1))["ll","Estimate"])
	PInc.2 <- exp(coef(summary(bm.2))["ll","Estimate"])

	#	increased odds - s.e.
	PInc.1.se <- exp(coef(summary(bm.1))["ll","Std. Error"])
	PInc.2.se <- exp(coef(summary(bm.2))["ll","Std. Error"])

	#	populate data.frame
	res.1[i,] <- c(run,x,y,z,1,ct,ct.ll,predict.ct.ll.1$fit,predict.ct.ll.1$se.fit,coef(summary(bm.1))["(Intercept)","Estimate"],bm.1.int.p,coef(summary(bm.1))["ll","Estimate"],bm.1.slope.p,PInc.1,PInc.1.se)
	res.2[i,] <- c(run,x,y,z,2,ct,ct.ll,predict.ct.ll.2$fit,predict.ct.ll.2$se.fit,coef(summary(bm.2))["(Intercept)","Estimate"],bm.2.int.p,coef(summary(bm.2))["ll","Estimate"],bm.2.slope.p,PInc.2,PInc.2.se)
}

res <- rbind(res.1,res.2)

# export results
write.table(res,file="ll.csv",row.names=F,append=T,sep=",",col.names=F)


