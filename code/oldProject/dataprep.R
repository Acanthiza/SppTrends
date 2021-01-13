###	FIX NAMES

#
# Clear workspace
#

	rm(list = ls())

#
# Import and setup data
#

	require(RODBC)
	con <- odbcConnectAccess("LRG/Rel/LRG_rel.mdb")
	tLUSpp <- sqlFetch(con,"tLUSpp")
	dat <- sqlFetch(con,"q00010")
	odbcCloseAll()
	rm(con)

	head(dat)
	head(tLUSpp)

	names(dat) <- c("poly","hex","com","year","records","sites","class")
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

###
### Summarise output with fixed names
###

#
# Overall species list
#

	spplist <- data.frame(spp=levels(dat$com))

#
# Overall poly list
#

	polylist <- data.frame(poly=levels(dat$poly))

#
# Maximum and minimum year for each species*landscape
#

	spppolymaxyear <- aggregate(year ~ poly + com + type, data = dat, max)
	spppolyminyear <- aggregate(year ~ poly + com + type, data = dat, min)

#
# Number of hexagons with a record for each species*landscape, overall ...
#

	spppolyhex <- aggregate(year ~ poly + com + type + hex, data = dat, length)

#
# Create summary dataframe including number of hexagons and max/min year
#

	datsum <- aggregate(hex ~ poly + com + type, data=spppolyhex, length)
	datsum <- cbind(datsum,maxyear=spppolymaxyear$year,minyear=spppolyminyear$year)

#
# Replace NA with 0
#

	datsum[is.na(datsum)] <- 0

	head(datsum)
	tail(datsum)

###
### Export data to .csv
###

#
# ensure directory for dataprep code table outputs
#

	dir.create(file.path("LRG","tbl","dataprep"))

#
# export spp list
#

	write.table(spplist,file=file.path("LRG","tbl","dataprep","spplist.csv"),row.names=F,append=F,sep=",",col.names=T)

#
# export spp*ls list
#
	
	write.table(datsum, file=file.path("LRG","tbl","dataprep","spppoly.csv"),row.names=F,append=F,sep=",",col.names=T)


#
# export poly list
#
	
	write.table(polylist, file=file.path("LRG","tbl","dataprep","polylist.csv"),row.names=F,append=F,sep=",",col.names=T)

#
# export results
#

	write.table(dat,file=file.path("LRG","tbl","dataprep","raw_withfixednames.csv"),row.names=F,append=F,sep=",",col.names=T)

#
# plot histogram of species records vs. year
#

	png(file=file.path("LRG","figs","dataprep","hist.png")
			, width=2000, height=2000, res=200
			, family = "Century Gothic"
			, bg = "transparent"
			)
	
	hist(bigdat$year
		, xlab = "Year"
		, main = ""
		)
	
	dev.off()


