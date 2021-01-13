#
# clear R workspace
#

	rm(list = ls())

#
# load required packages
#

	library(rgdal)			# import spatial data (readOGR)
	library(adehabitatHR)		# run locoh analyses

#
# get data from personal geodatabase
#
	#TEST - includes 'Sub'
	spdat <- readOGR(file.path("LRG","gdb","LRG_Geo.mdb"),"STHexSub")

	names(spdat) <- c("oid","lat","lon","com","year")

#
# OR load workspace
#

	#load(file.path("LRG","R","locoh_shp.RData"))

#
# fix spp names
#

	spdat$com <- gsub("'Adelaide Rosella'","Crimson Rosella",spdat$com,fixed=T)
	spdat$com <- gsub("Australian Ringneck, (Ring-necked Parrot)","Australian Ringneck",spdat$com,fixed=T)
	spdat$com <- gsub("Chestnut-rumped Heathwren (ML Ranges ssp)","Chestnut-rumped Heathwren",spdat$com,fixed=T)
	spdat$com <- gsub("Chestnut-rumped Heathwren (South East ssp)","Chestnut-rumped Heathwren",spdat$com,fixed=T)
	spdat$com <- gsub("Chestnut Quail-thrush (eastern ssp)","Chestnut Quailthrush",spdat$com,fixed=T)
	spdat$com <- gsub("Clamorous Reedwarbler","Australian Reed Warbler",spdat$com,fixed=T)
	spdat$com <- gsub("Glossy Black-Cockatoo (Kangaroo Island ssp)","Glossy Black-Cockatoo",spdat$com,fixed=T)
	spdat$com <- gsub("Hooded Robin (South East ssp)","Hooded Robin",spdat$com,fixed=T)
	spdat$com <- gsub("Pacific Black Duck/Mallard Hybrid","Pacific Black Duck",spdat$com,fixed=T)
	spdat$com <- gsub("Port Lincoln Parrot","Australian Ringneck",spdat$com,fixed=T)
	spdat$com <- gsub("Red-tailed Black Cockatoo (south-east subspecies)","Red-tailed Black Cockatoo",spdat$com,fixed=T)
	spdat$com <- gsub("Slender-billed Thornbill (western ssp)","Slender-billed Thornbill",spdat$com,fixed=T)
	spdat$com <- gsub("Southern Emu-wren (Mt Lofty Ranges ssp)","Southern Emuwren",spdat$com,fixed=T)
	spdat$com <- gsub("Southern Emu-wren","Southern Emuwren",spdat$com,fixed=T)
	spdat$com <- gsub("Southern Emu-wren (South East ssp)","Southern Emuwren",spdat$com,fixed=T)
	spdat$com <- gsub("Southern Emuwren (Eyre Peninsula ssp)","Southern Emuwren",spdat$com,fixed=T)
	spdat$com <- gsub("Southern Emuwren (Kangaroo Island ssp)","Southern Emuwren",spdat$com,fixed=T)
	spdat$com <- gsub("Southern Emuwren (South East ssp)","Southern Emuwren",spdat$com,fixed=T)
	spdat$com <- gsub("Spotted Quail-thrush (Mount Lofty Ranges ssp)","Spotted Quailthrush",spdat$com,fixed=T)
	spdat$com <- gsub("Spotted Quail-thrush","Spotted Quailthrush",spdat$com,fixed=T)
	spdat$com <- gsub("Spur-winged Plover","Masked Lapwing",spdat$com,fixed=T)
	spdat$com <- gsub("Western Whipbird (Eastern subspecies)","Western Whipbird",spdat$com,fixed=T)
	spdat$com <- gsub("Western Whipbird (Kangaroo Island ssp)","Western Whipbird",spdat$com,fixed=T)
	spdat$com <- gsub("Yellow Rosella","Crimson Rosella",spdat$com,fixed=T)
	spdat$com <- gsub("Yellow-tailed Pardalote","Spotted Pardalote",spdat$com,fixed=T)
	spdat$com <- gsub("Yellow-throated/Black-eared Miner Cross","Yellow-throated Miner",spdat$com,fixed=T)
	spdat$com <- gsub("Yellow-vented Bluebonnet","Bluebonnet",spdat$com,fixed=T)
	spdat$com <- gsub("Slender-billed Thornbill (St Vincent Gulf ssp)","Slender-billed Thornbill",spdat$com,fixed=T)
	spdat$com <- gsub("Brown Hawk (Brown Falcon)","Brown Falcon",spdat$com,fixed=T)
	spdat$com <- gsub("Naretha Bluebonnet","Bluebonnet",spdat$com,fixed=T)

#
# clean up data set
#

	remove.incomplete.cases <- which(complete.cases(spdat@data) == FALSE)
	spdat <- spdat[-remove.incomplete.cases,]

#
# split into list by species
#

	spdat.split <- split(spdat, spdat$com)

#
# remove spp from list with < 20 records
#

	spp <- aggregate(oid ~ com, data = spdat, FUN = length)
	spp.remove <- which(spp[,2] < 20)
	spdat.split <- spdat.split[-spp.remove]

#
# remove points closer to each other than 500m, then reassemble spatial data.frame as spdat.filt
#
	
	spdat.filt <- spdat[0,]

	for (j in 1:length(spdat.split)){
		out <- spdat.split[[j]][1,]
		for(i in 2:nrow(spdat.split[[j]])){
			pt <- coordinates(spdat.split[[j]][i,])
			pts <-  coordinates(spdat.split[[j]][-i,])
			dists <- spDistsN1(pts, pt, longlat = F)
			exceed <- any(dists < 500)
			if (exceed==FALSE){
				out <- rbind(out, spdat.split[[j]][i,])
				}
			}
		spdat.filt <- rbind(out, spdat.filt)
		}

#
# 
#

	use.oid <- data.frame(oid=spdat.filt$oid)

	spdat.sub <- spdat[spdat$oid  use.oid,]
	
	areaKcurve <- LoCoH.k.area(spdat.filt.sub[,4]
		, k = c(5:15)
		, percent=100
		)



