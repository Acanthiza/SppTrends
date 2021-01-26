
  

  #-------AOI---------
  aoi <- make_aoi(polys
                  ,filterPolys = aoiName
                  ,filterPolysCol = "LSA"
                  ,polyBuffer = polyBuffer
                  )
  
  sitesAOI <- sitesAll %>%
    st_join(aoi
            , join = st_intersects
            ) %>%
    dplyr::filter(!is.na(Include))

  taxaAllAOI <- taxaAll %>%
    dplyr::inner_join(sitesAOI %>% st_set_geometry(NULL) %>% dplyr::select(LATITUDE,LONGITUDE))
  
  
  #-------Taxonomy-------
  
  taxaAllAOITax <- taxaAllAOI %>%
    dplyr::left_join(luGBIF[,c("id","Taxa","Rank","Class")], by = c("SPECIES" = "id")) %>%
    dplyr::filter(Rank > "Genus") %>%
    dplyr::count(Class,Taxa,LATITUDE,LONGITUDE,year,name="siteRecords")
  
  #-------Add geo context-------
  
  patchGeo <- taxaAllAOITax %>%
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
    dplyr::rename(IBRASub = IBRA_SUB_N) %>%
    dplyr::distinct(LATITUDE,LONGITUDE,cell,IBRASub)
  
  taxaAllAOITaxGeo <- taxaAllAOITax %>%
    dplyr::left_join(patchGeoContext) %>%
    dplyr::distinct(Class,Taxa,year,IBRASub,cell) %>%
    dplyr::mutate(list = paste0(year,"-",Class,"-",cell))
  
  
  #------Filter--------
  
  listLength <- taxaAllAOITaxGeo %>%
    dplyr::count(list,name="listLength") %>%
    dplyr::filter(listLength > minRRListLength)
  
  lists <- taxaAllAOITaxGeo %>%
    dplyr::group_by(Class,year,IBRASub) %>%
    dplyr::summarise(lists = n_distinct(cell)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(lists > 3)
  
  yearsIBRASub <- taxaAllAOITaxGeo %>%
    dplyr::group_by(IBRASub) %>%
    dplyr::summarise(yearsIBRASub = n_distinct(year)) %>%
    dplyr::ungroup()
  
  yearsIBRASubTaxa <- taxaAllAOITaxGeo %>% 
    dplyr::group_by(IBRASub,Taxa) %>%
    dplyr::summarise(yearsIBRASubTaxa = n_distinct(year)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(yearsIBRASubTaxa > minTaxaOccurence)
    
  taxaAllAOITaxGeoFilter <- taxaAllAOITaxGeo %>%
    dplyr::inner_join(listLength) %>%
    dplyr::inner_join(lists) %>%
    dplyr::inner_join(yearsIBRASub) %>%
    dplyr::inner_join(yearsIBRASubTaxa)
    
    
  
  