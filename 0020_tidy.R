
  timer$start("tidy")
  
  alwaysGroup <- c("year","yday","month","geo1","geo2","cell","site") #in order
  

  #------Dates-------
  
  taxaAllDates <- taxaAll %>%
    dplyr::filter(year >= minYear) %>%
    dplyr::mutate(quart = cut(month
                              , breaks = c(0,4,6,9,12)
                              , labels = c("q1","q2","q3","q4")
                              )
                  )
  

  #-------AOI---------
  
  sitesAll <- taxaAllDates %>%
    dplyr::distinct(LATITUDE,LONGITUDE) %>%
    st_as_sf(coords = c("LONGITUDE","LATITUDE"), crs = 4326, remove = FALSE) %>%
    st_transform(crs = st_crs(polys))
  
  aoi <- make_aoi(polygons = polys
                  ,filterPolys = polyMask
                  ,filterPolysCol = "IBRA_REG_N"
                  ,polyBuffer = polyBuf
                  )
  
  sitesAOI <- sitesAll %>%
    sf::st_join(aoi
                , left = FALSE
                )

  taxaAllDatesAOI <- taxaAllDates %>%
    dplyr::inner_join(sitesAOI %>% st_set_geometry(NULL) %>% dplyr::select(LATITUDE,LONGITUDE))
  
  
  #-------Taxonomy-------
  
  taxaAllDatesAOITax <- taxaAllDatesAOI %>%
    dplyr::left_join(luGBIF[,c("id","Taxa","Rank")], by = c("SPECIES" = "id")) %>%
    dplyr::filter(Rank > "Genus") %>%
    dplyr::count(Taxa,LATITUDE,LONGITUDE,date,year,month,yday,quart,yearmon,name="siteRecords") %>%
    dplyr::left_join(luTax) %>%
    dplyr::filter(!!ensym(taxGroup) %in% testTaxGroups)
  
  
  #-------Add geo context-------
  
  siteGeo <- taxaAllDatesAOITax %>%
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
    dplyr::distinct(LATITUDE,LONGITUDE,site,cell,geo1,geo2) %>%
    dplyr::filter(geo1 %in% polyMask)
  
  taxaAllDatesAOITaxGeo <- taxaAllDatesAOITax %>%
    dplyr::inner_join(siteGeoContext) %>%
    add_list_length(groups = c(taxGroup,alwaysGroup))
  
  
  #-------Cooccuring Taxa-----------
  
  contextCooccur <- taxaAllDatesAOITaxGeo %>%
    dplyr::filter(listLength > 1) %>%
    tidyr::nest(data = -c(!!ensym(taxGroup),geo1,geo2)) %>%
    dplyr::mutate(nLists = map_dbl(data,~n_distinct(.$list))
                  , nTaxa = map_dbl(data,~n_distinct(.$Taxa))
                  ) %>%
    # more than 10*min(possibleGroups) lists within a context
    dplyr::filter(nLists > 10*min(possibleGroups)) %>%
    # at least 2*min(possibleGroups) distinct Taxa
    dplyr::filter(nTaxa > 2*min(possibleGroups)) %>%
    # at least three taxa with more than three records
    dplyr::filter(map_dbl(data,~sum(table(.$Taxa) < 3)) > 3) %>%

    #dplyr::sample_n(5) %>% # TESTING

    dplyr::mutate(cooccur = future_map(data,make_cooccur))
  
  sppCooccur <- contextCooccur %>%
    dplyr::mutate(pair = map(cooccur,"results")) %>%
    tidyr::unnest(cols = c(pair)) %>%
    dplyr::select(where(negate(is.list))) %>%
    dplyr::filter(p_gt < 0.05) %>%
    dplyr::rename(Taxa = sp1_name)
    
  
  #------datTidy--------
  
  datTidy <- taxaAllDatesAOITaxGeo %>%
    dplyr::mutate(cell = as.factor(cell)
                  , site = as.factor(site)
                  ) %>%
    dplyr::left_join(luTax[,c("Taxa","Common")])
  
  
  timer$stop("tidy", comment = paste0("tidy took records from ",nrow(taxaAll)," to ",nrow(datTidy)))
  