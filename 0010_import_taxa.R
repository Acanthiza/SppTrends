
  timer$start("import")
  
  #--------BDBSA---------
  
  source(path("code","bdbsa.r"))
  
  taxaBDBSA <- rawBDBSA %>%
    dplyr::mutate(year = year(VISITDATE)
                  , month = month(VISITDATE)
                  , yearmon = as.numeric(paste0(year,sprintf("%02d",month)))
                  ) %>%
    dplyr::add_count(SPECIES, name = "records") %>%
    dplyr::filter(!is.na(LATITUDE)
                  , !is.na(LONGITUDE)
                  , !is.na(SPECIES)
                  , !is.na(yearmon)
                  , !grepl("in-active|sfossil",METHODDESC)
                  , !is.na(NSXCODE)
                  , !NUMOBSERVED %in% nonRecords
                  , records > 3
                  ) %>%
    dplyr::rename(CommonName = COMNAME1) %>%
    dplyr::left_join(luRel) %>%
    dplyr::select(any_of(collectFields)) %>%
    dplyr::mutate(source = "BDBSA")
  
  
  #---------GBIF---------
  
  gbifFile <- path("out","gbif","taxaGBIF.feather")
  
  if(!file.exists(gbifFile)) {
    
    source(path("code","gbif.r"))
  
    taxaGBIF <- rawGBIF %>%
      as_tibble() %>%
      dplyr::mutate(SPECIES = str_extract(scientificName,"[[:alpha:]]+\\s[[:alpha:]]+")
                    , year = as.numeric(substr(eventDate,1,4))
                    , month = as.numeric(substr(eventDate,6,7))
                    , yearmon = as.numeric(paste0(year,substr(eventDate,6,7)))
                    ) %>%
      dplyr::select(LATITUDE = decimalLatitude
                    , LONGITUDE = decimalLongitude
                    , SPECIES
                    , CommonName = vernacularName
                    #, METHODDESC = ?
                    , NUMOBSERVED = individualCount
                    , year
                    , month
                    , yearmon
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
    
    write_feather(taxaGBIF,gbifFile)
    
  }
  
  taxaGBIF <- read_feather(gbifFile)
  
  
  #-------Combine-------
  taxaAll <- ls(pattern = "^taxa[[:upper:]]+") %>%
    enframe(name = NULL, value = "objects") %>%
    dplyr::filter(!grepl("All",objects)) %>%
    dplyr::mutate(data = map(objects,get)) %>%
    tidyr::unnest(cols = c(data)) %>%
    dplyr::mutate(SPECIES = gsub("\\s"," ",SPECIES))
  
  
  
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
                                  , ordered = TRUE
                                  )
                  ) %>%
    dplyr::mutate(Taxa = if_else(originalName == "Pandion haliaetus","Pandion cristatus",Taxa)
                  , Common = if_else(grepl("^Osprey$",Common),"Eastern Osprey",Common)
                  )
  
  luTax <- luGBIF %>%
    dplyr::distinct(Taxa,Common,Kingdom,Phylum,Class,Order,Family,Genus)
  
  luInd <- taxaBDBSA %>%
    dplyr::left_join(luGBIF, by = c("SPECIES" = "id")) %>%
    dplyr::filter(!is.na(Taxa)) %>%
    dplyr::count(Taxa,ISINDIGENOUSFLAG) %>%
    dplyr::mutate(ISINDIGENOUS = if_else(is.na(ISINDIGENOUSFLAG),"Y","N")) %>%
    dplyr::group_by(Taxa) %>%
    dplyr::filter(n == max(n, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::select(Taxa,ISINDIGENOUS)
  
  
  timer$stop("import", comment = paste0("Get new data = ",getNewData))
  