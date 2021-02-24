
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
  
  
  #-------taxGroup---------
  
  taxGroups <- luTax %>%
    dplyr::filter(Taxa %in% tests) %>%
    dplyr::pull(!!ensym(taxGroup)) %>%
    unique()
  
  taxaAllDatesAOITaxTaxGroup <- taxaAllDatesAOITax %>%
    dplyr::filter(!!ensym(taxGroup) %in% taxGroups)
  
  
  #-------Add geo context-------
  
  siteGeo <- taxaAllDatesAOITaxTaxGroup %>%
    dplyr::distinct(LONGITUDE,LATITUDE) %>%
    st_as_sf(coords = c("LONGITUDE","LATITUDE"), crs = 4326, remove = FALSE) %>%
    st_transform(crs = crs(polys))
  
  siteExtent <- extent(siteGeo)
  
  innerRaster <- raster(siteExtent
                      , res = c(innerGrid,innerGrid)
                      , crs = crs(polys)
                      )
  
  outerRaster <- raster(siteExtent
                      , res = c(outerGrid,outerGrid)
                      , crs = crs(polys)
                      )
  
  siteGeoContext <- siteGeo %>%
    dplyr::mutate(site = cellFromXY(innerRaster,as_Spatial(siteGeo$geometry))
                  , cell = cellFromXY(outerRaster,as_Spatial(siteGeo$geometry))
                  ) %>%
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
                  , !is.na(site)
                  ) %>%
    dplyr::distinct(LATITUDE,LONGITUDE,site,cell,geo1,geo2)
  
  
  #------dat--------
  datTidy <- taxaAllDatesAOITaxTaxGroup %>%
    dplyr::inner_join(siteGeoContext) %>%
    dplyr::distinct(!!ensym(taxGroup),Taxa,year,month,yearmon,geo1,geo2,cell,site) %>%
    dplyr::mutate(cell = as.factor(cell)
                  , site = as.factor(site)
                  )
  
  
  timer$stop("tidy", comment = paste0("tidy took records from ",nrow(taxaAll)," to ",nrow(datTidy)))
  