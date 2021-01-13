#
# clear R workspace
#

	rm(list = ls())

#
# load required packages
#

	library(adehabitatHR)
	library(RODBC)
	library(maptools)			# import spatial data
	#library(RgoogleMaps)
	#library(ggmap)	

#
# set up data
#

	con <- odbcConnectAccess("LRG/Rel/LRG_rel.mdb")
	tLUSpp <- sqlFetch(con,"tLUSpp")
	dat <- sqlFetch(con,"q00020")
	odbcCloseAll()
	rm(con)

	head(dat)
	head(tLUSpp)

	names(dat) <- c("poly","long","lat","hex","com","year","records","sites","class")
	names(tLUSpp) <-c("type","com","n","class","native")
	tLUSpp$native <- as.factor(tLUSpp$native)

	#dat$com <- gsub("","",dat$com)

	dat$com <- gsub("'Adelaide Rosella'","Crimson Rosella",dat$com,fixed=T)
	dat$com <- gsub("Australian Ringneck, (Ring-necked Parrot)","Australian Ringneck",dat$com,fixed=T)
	dat$com <- gsub("Chestnut-rumped Heathwren (ML Ranges ssp)","Chestnut-rumped Heathwren",dat$com,fixed=T)
	dat$com <- gsub("Chestnut-rumped Heathwren (South East ssp)","Chestnut-rumped Heathwren",dat$com,fixed=T)
	dat$com <- gsub("Chestnut Quail-thrush (eastern ssp)","Chestnut Quailthrush",dat$com,fixed=T)
	dat$com <- gsub("Clamorous Reedwarbler","Australian Reed Warbler",dat$com,fixed=T)
	dat$com <- gsub("Glossy Black-Cockatoo (Kangaroo Island ssp)","Glossy Black-Cockatoo",dat$com,fixed=T)
	dat$com <- gsub("Hooded Robin (South East ssp)","Hooded Robin",dat$com,fixed=T)
	dat$com <- gsub("Pacific Black Duck/Mallard Hybrid","Pacific Black Duck",dat$com,fixed=T)
	dat$com <- gsub("Port Lincoln Parrot","Australian Ringneck",dat$com,fixed=T)
	dat$com <- gsub("Red-tailed Black Cockatoo (south-east subspecies)","Red-tailed Black Cockatoo",dat$com,fixed=T)
	dat$com <- gsub("Slender-billed Thornbill (western ssp)","Slender-billed Thornbill",dat$com,fixed=T)
	dat$com <- gsub("Southern Emu-wren (Mt Lofty Ranges ssp)","Southern Emuwren",dat$com,fixed=T)
	dat$com <- gsub("Southern Emu-wren","Southern Emuwren",dat$com,fixed=T)
	dat$com <- gsub("Southern Emu-wren (South East ssp)","Southern Emuwren",dat$com,fixed=T)
	dat$com <- gsub("Southern Emuwren (Eyre Peninsula ssp)","Southern Emuwren",dat$com,fixed=T)
	dat$com <- gsub("Southern Emuwren (Kangaroo Island ssp)","Southern Emuwren",dat$com,fixed=T)
	dat$com <- gsub("Southern Emuwren (South East ssp)","Southern Emuwren",dat$com,fixed=T)
	dat$com <- gsub("Spotted Quail-thrush (Mount Lofty Ranges ssp)","Spotted Quailthrush",dat$com,fixed=T)
	dat$com <- gsub("Spotted Quail-thrush","Spotted Quailthrush",dat$com,fixed=T)
	dat$com <- gsub("Spur-winged Plover","Masked Lapwing",dat$com,fixed=T)
	dat$com <- gsub("Western Whipbird (Eastern subspecies)","Western Whipbird",dat$com,fixed=T)
	dat$com <- gsub("Western Whipbird (Kangaroo Island ssp)","Western Whipbird",dat$com,fixed=T)
	dat$com <- gsub("Yellow Rosella","Crimson Rosella",dat$com,fixed=T)
	dat$com <- gsub("Yellow-tailed Pardalote","Spotted Pardalote",dat$com,fixed=T)
	dat$com <- gsub("Yellow-throated/Black-eared Miner Cross","Yellow-throated Miner",dat$com,fixed=T)
	dat$com <- gsub("Yellow-vented Bluebonnet","Bluebonnet",dat$com,fixed=T)
	dat$com <- gsub("Slender-billed Thornbill (St Vincent Gulf ssp)","Slender-billed Thornbill",dat$com,fixed=T)
	dat$com <- gsub("Brown Hawk (Brown Falcon)","Brown Falcon",dat$com,fixed=T)
	dat$com <- gsub("Naretha Bluebonnet","Bluebonnet",dat$com,fixed=T)

	dat <- merge(dat,tLUSpp[,c(1,2,5)])

	dat$native <- as.factor(dat$native)

#
# dump fixed names to file (for easy reload, the access link takes a while)
#

	write.table(dat
		, file=file.path("LRG","tbl","dataprep","locoh_data.csv")
		, row.names=F
		, append=F
		, sep=","
		, col.names=T
		)

#
# setup data
#

	#dat <- read.csv(file.path("LRG","tbl","dataprep","locoh_data.csv"))

	dat$long <- round(dat$long,3)
	dat$lat <- round(dat$lat,3)
	dat <- aggregate(type ~ com + long + lat, data = dat, FUN = length)
	
	dat.split <- split(dat, dat$com)
	spp <- aggregate(type ~ com, data = dat, FUN = length)
	spp.remove <- which(spp$type < 10)
	dat.split <- dat.split[-spp.remove]

	dat.filter <- dat[0,]

	for (j in 1:length(dat.split)){
		out <- dat.split[[j]][1,]
		for(i in 2:nrow(dat.split[[j]])){
			pts <- as.matrix(out[,2:3])
			pt <-  as.numeric(dat.split[[j]][i,2:3])
			dists <- spDistsN1(pts, pt, longlat = T)
			exceed <- any(dists < 0.5)
			if (exceed==FALSE){
				out <- rbind(out, dat.split[[j]][i,])
				}
			}
		dat.filter <- rbind(out, dat.filter)
		}
		

	dat.sub <- dat[dat$com == "Chestnut Quailthrush",]

	dat.crs <- CRS("+proj=longlat +datum=WGS84")

# see http://www.spatialreference.org/ref/epsg/3107/

	points <- dat.sub[,2:3]
	points <- SpatialPoints(points
			, proj4string = dat.crs
		)	

	areaKcurve <- LoCoH.k.area(points
		, k = c(8:12)
		, percent=90
		#, unin = "m"
		#, unout = "ha"
		, duplicates = "remove"
		#, amount = NULL
		)

