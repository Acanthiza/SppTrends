
###	Prepare mapping query
## working directory
setwd("C:/Workspace_NW/LRG/shp/")

##Import and setup data
rm(list=ls())
require(RODBC)
con <- odbcConnectAccess("C:/Workspace_NW/LRG/rel/LRG_Prep_rel.mdb")
dat <- sqlFetch(con,"q01470")
odbcCloseAll()
rm(con)

head(dat)

# export results
dat$MapID <- seq(1:nrow(dat))

require(shapefiles)

dat.xy <- data.frame(MapID=dat$MapID,E=dat$e,N=dat$n)
dat.tab <- data.frame(MapID=dat$MapID,Poly=dat$poly,Name=dat$name,Decl=dat$Decl,Dist=dat$Dist,Exti=dat$Exti,Year=dat$year,records=dat$records)
dat.shp <- convert.to.shapefile(dat.xy, dat.tab, "MapID", 1)
write.shapefile(dat.shp, "Map_LRG", arcgis=T)


