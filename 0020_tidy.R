
  timer$start("tidy")

  #------Dates-------
  taxaAllDates <- taxaAll %>%
    dplyr::filter(year >= minYear)

  #-------AOI---------
  sitesAll <- taxaAllDates %>%
    dplyr::distinct(LATITUDE,LONGITUDE) %>%
    st_as_sf(coords = c("LONGITUDE","LATITUDE"), crs = 4326, remove = FALSE) %>%
    st_transform(crs = st_crs(polys))
  
  aoi <- make_aoi(polys
                  ,filterPolys = aoiName
                  ,filterPolysCol = "IBRA_REG_C"
                  ,polyBuffer = polyBuffer
                  )
  
  sitesAOI <- sitesAll %>%
    st_join(aoi
            , join = st_intersects
            ) %>%
    dplyr::filter(!is.na(Include))

  taxaAllDatesAOI <- taxaAllDates %>%
    dplyr::inner_join(sitesAOI %>% st_set_geometry(NULL) %>% dplyr::select(LATITUDE,LONGITUDE))
  
  
  #-------Taxonomy-------
  
  taxaAllDatesAOITax <- taxaAllDatesAOI %>%
    dplyr::left_join(luGBIF[,c("id","Taxa","Rank")], by = c("SPECIES" = "id")) %>%
    dplyr::filter(Rank > "Genus") %>%
    dplyr::count(Taxa,LATITUDE,LONGITUDE,year,month,yearmon,name="siteRecords") %>%
    dplyr::left_join(luTax)
  
  #-------Add geo context-------
  
  patchGeo <- taxaAllDatesAOITax %>%
    dplyr::distinct(LONGITUDE,LATITUDE) %>%
    st_as_sf(coords = c("LONGITUDE","LATITUDE"), crs = 4326, remove = FALSE) %>%
    st_transform(crs = crs(polys))
  
  patchExtent <- extent(patchGeo)
  
  patchRaster <- raster(patchExtent
                       , res = c(gridSize,gridSize)
                       , crs = crs(polys)
                       )
  
  patchGeoContext <- patchGeo %>%
    dplyr::mutate(cell = cellFromXY(patchRaster,as_Spatial(patchGeo$geometry))) %>%
    st_join(ibraSub
            , join = st_intersects
            ) %>%
    st_set_geometry(NULL) %>%
    as_tibble() %>%
    dplyr::rename(geo1 = IBRA_REG_N
                  , geo2 = IBRA_SUB_N
                  ) %>%
    dplyr::filter(!is.na(geo1)
                  , !is.na(geo2)
                  , !is.na(cell)
                  ) %>%
    dplyr::distinct(LATITUDE,LONGITUDE,cell,geo1,geo2)
  
  
  #------dat--------
  dat <- taxaAllDatesAOITax %>%
    dplyr::inner_join(patchGeoContext) %>%
    dplyr::distinct(!!ensym(taxGroup),Taxa,year,month,yearmon,geo1,geo2,cell) %>%
    dplyr::mutate(cell = as.factor(cell))
  
  
  timer$stop("tidy", comment = paste0("tidy took records from ",nrow(taxaAll)," to ",nrow(dat)))