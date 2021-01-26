
  #--------BDBSA---------
  
  source(path("code","bdbsa.r"))
  
  taxaBDBSA <- rawBDBSA %>%
    dplyr::mutate(year = year(VISITDATE)) %>%
    dplyr::add_count(SPECIES, name = "records") %>%
    dplyr::filter(!is.na(LATITUDE)
                  , !is.na(LONGITUDE)
                  , !is.na(SPECIES)
                  , !is.na(year)
                  #, !grepl("in-active|sfossil",METHODDESC)
                  , !is.na(NSXCODE)
                  , !NUMOBSERVED %in% nonRecords
                  , records > 3
                  ) %>%
    dplyr::rename(CommonName = COMNAME1) %>%
    dplyr::left_join(luRel) %>%
    dplyr::select(any_of(collectFields)) %>%
    dplyr::mutate(source = "BDBSA")
  
  
  #---------GBIF---------
  
  source(path("code","gbif.r"))
  
  taxaGBIF <- rawGBIF %>%
    as_tibble() %>%
    dplyr::mutate(SPECIES = str_extract(scientificName,"[[:alpha:]]+\\s[[:alpha:]]+")
                  , year = as.numeric(substr(year,1,4))
                  ) %>%
    dplyr::select(LATITUDE = decimalLatitude
                  , LONGITUDE = decimalLongitude
                  , SPECIES
                  , CommonName = vernacularName
                  #, METHODDESC = ?
                  , NUMOBSERVED = individualCount
                  , year
                  , maxDist = coordinateUncertaintyInMeters
                  ) %>%
    dplyr::add_count(SPECIES, name = "records") %>%
    dplyr::filter(!is.na(LATITUDE)
                  , !is.na(LONGITUDE)
                  , !is.na(SPECIES)
                  , !is.na(year)
                  , records > 3
                  ) %>%
    dplyr::mutate(source = "GBIF"
                  , NUMOBSERVED = as.character(NUMOBSERVED)
                  , maxDist = gsub('>|""',NA,maxDist)
                  , maxDist = as.numeric(maxDist)
                  )
  
  
  #-------Combine-------
  taxaAll <- ls(pattern = "^taxa[[:upper:]]+") %>%
    enframe(name = NULL, value = "objects") %>%
    dplyr::filter(!grepl("All",objects)) %>%
    dplyr::mutate(data = map(objects,get)) %>%
    tidyr::unnest(cols = c(data)) %>%
    dplyr::mutate(SPECIES = gsub("\\s"," ",SPECIES))
  
  sitesAll <- taxaAll %>%
    dplyr::distinct(LATITUDE,LONGITUDE) %>%
    st_as_sf(coords = c("LONGITUDE","LATITUDE"), crs = 4326, remove = FALSE) %>%
    st_transform(crs = st_crs(polys))
  
  #-------Taxonomy------
  taxaAll %>%
    dplyr::distinct(SPECIES) %>%
    #dplyr::sample_n(30) %>% # FOR TESTING ONLY
    gbif_tax(1,1,path("data","luGBIF.feather"),"Animalia")
  
  luGBIF <- read_feather(path("data","luGBIF.feather")) %>%
    dplyr::mutate(TaxaGBIF = Taxa
                  , Rank = str_to_sentence(Rank)
                  , Rank = gsub("Division","Phylum",Rank)
                  , Rank = fct_expand(Rank,"Class","Order")
                  , Rank = factor(Rank
                                  , levels = c("Kingdom","Phylum","Family","Genus","Species","Subspecies","Variety","Form")
                                  , ordered = TRUE)
                  )
  
  luTax <- luGBIF %>%
    dplyr::distinct(Taxa,Class)
  
  luInd <- taxaBDBSA %>%
    dplyr::left_join(luGBIF, by = c("SPECIES" = "id")) %>%
    dplyr::filter(!is.na(Taxa)) %>%
    dplyr::count(Taxa,ISINDIGENOUSFLAG) %>%
    dplyr::mutate(ISINDIGENOUS = if_else(is.na(ISINDIGENOUSFLAG),"Y","N")) %>%
    dplyr::group_by(Taxa) %>%
    dplyr::filter(n == max(n, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::select(Taxa,ISINDIGENOUS)
  