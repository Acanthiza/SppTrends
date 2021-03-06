
  timer$start("import")
  
  #--------BDBSA---------
  
  source(path("code","bdbsa.r"))
  
  taxaBDBSA <- rawBDBSA %>%
    dplyr::mutate(date = VISITDATE
                  , year = year(VISITDATE)
                  , month = month(VISITDATE)
                  , yday = yday(VISITDATE)
                  , yearmon = as.numeric(paste0(year,sprintf("%02d",month)))
                  ) %>%
    dplyr::add_count(SPECIES, name = "records") %>%
    dplyr::filter(!is.na(LATITUDE)
                  , !is.na(LONGITUDE)
                  , !is.na(date)
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
      dplyr::mutate(date = ymd(gsub("T00:00:00","",eventDate))
                    , SPECIES = scientificName
                    , year = as.numeric(substr(eventDate,1,4))
                    , month = as.numeric(substr(eventDate,6,7))
                    , yday = yday(date)
                    , yearmon = as.numeric(paste0(year,substr(eventDate,6,7)))
                    ) %>%
      dplyr::distinct(LATITUDE = decimalLatitude
                    , LONGITUDE = decimalLongitude
                    , date
                    , year
                    , month
                    , yday
                    , yearmon
                    , SPECIES
                    , CommonName = vernacularName
                    #, METHODDESC = ?
                    , NUMOBSERVED = individualCount
                    , maxDist = coordinateUncertaintyInMeters
                    ) %>%
      dplyr::add_count(SPECIES, name = "records") %>%
      dplyr::filter(!is.na(LATITUDE)
                    , !is.na(LONGITUDE)
                    , !is.na(SPECIES)
                    , !is.na(date)
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
  
  taxaAll <- ls(pattern = "^taxa[[:upper:]]+$") %>%
    enframe(name = NULL, value = "objects") %>%
    dplyr::filter(!grepl("All",objects)) %>%
    dplyr::mutate(data = map(objects,get)) %>%
    tidyr::unnest(cols = c(data)) %>%
    dplyr::mutate(SPECIES = gsub("\\s"," ",SPECIES))
  
  
  #-------Taxonomy------
  
  taxaAll %>%
    dplyr::distinct(SPECIES) %>%
    dplyr::arrange(SPECIES) %>%
    #dplyr::sample_n(30) %>% # FOR TESTING ONLY
    gbif_tax(1
             ,path("data","luGBIF.feather")
             ,"Animalia"
             ,getCommon = TRUE
             )
  
  luGBIF <- read_feather(path("data","luGBIF.feather")) %>%
    dplyr::mutate(TaxaGBIF = Taxa
                  , Rank = rank
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
    dplyr::distinct(Taxa,kingdom,phylum,class,order,family,genus,species,Common) %>%
    na.omit() %>%
    dplyr::arrange(Taxa)
  
  luInd <- taxaBDBSA %>%
    dplyr::left_join(luGBIF, by = c("SPECIES" = "originalName")) %>%
    dplyr::filter(!is.na(Taxa)) %>%
    dplyr::count(Taxa,ISINDIGENOUSFLAG) %>%
    dplyr::mutate(ISINDIGENOUS = if_else(is.na(ISINDIGENOUSFLAG),"Y","N")) %>%
    dplyr::group_by(Taxa) %>%
    dplyr::filter(n == max(n, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::select(Taxa,ISINDIGENOUS)
  
  
  timer$stop("import", comment = paste0("Get new data = ",getNewData))
  