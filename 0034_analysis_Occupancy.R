

  timer$start("occupancy")
  
  #-------datForOcc---------
  
  # Sampling unit for occupancy is based on a cell within a yearmon. Thus, success is presence within a cell in a yearmon
  
  allScales <- c("geo1","geo2","cell") #in order
  analysisScales <- allScales#[-length(allScales)]
  
  trials <- datTidy %>%
    dplyr::distinct(year,quart,!!ensym(taxGroup),across(any_of(allScales))) %>%
    dplyr::group_by(year,quart,!!ensym(taxGroup),across(any_of(analysisScales))) %>%
    dplyr::summarise(trials = n()) %>%
    dplyr::ungroup()
  
  success <- datTidy %>%
    dplyr::distinct(Taxa,year,quart,!!ensym(taxGroup),across(any_of(allScales))) %>%
    dplyr::group_by(Taxa,year,quart,!!ensym(taxGroup),across(any_of(analysisScales))) %>%
    dplyr::summarise(success = n()) %>%
    dplyr::ungroup()
  
  datFiltered <- trials %>%
    dplyr::left_join(success) %>%
    dplyr::mutate(list = paste0(year,"-",quart,"-",!!ensym(taxGroup),"-",get(analysisScales[length(analysisScales)]))) %>%
    dplyr::add_count(list, name = "listLength") %>%
    filter_taxa_data(minListLengthThresh = 0
                     , maxListLengthOccurenceThresh = 0
                     , minlistOccurenceThresh = 1
                     , minYearsThresh = 5
                     , minListLengthsThresh = 0
                     )
  
  taxaGeo <- datFiltered %>%
    dplyr::distinct(!!ensym(taxGroup)
                    ,Taxa
                    ,across(any_of(analysisScales))
                    )
  
  #-------datForRR---------
  
  dat <- datFiltered %>%
    dplyr::inner_join(taxaGeo) %>%
    tidyr::pivot_wider(names_from = c("Taxa","quart"), values_from = "success", values_fill = 0) %>%
    tidyr::pivot_longer(cols = contains(unique(taxaGeo$Taxa))
                        , names_to = c("Taxa","quart")
                        , names_sep = "_"
                        , values_to = "success"
                        ) %>%
    dplyr::inner_join(taxaGeo) %>%
    purrr::when(testing ~ (.) %>% dplyr::filter(Taxa %in% tests)
                , !testing ~ (.)
                ) %>%
    tidyr::nest(data = c(any_of(analysisScales),year,quart,list,listLength,success,trials)) %>%
    dplyr::left_join(luTax %>%
                       dplyr::select(Taxa,Common)
                     )
  
  #------Model function--------
  
  occ <- function(Taxa,data) {
    
    print(paste0(Taxa))
    
    outFile <- fs::path(outDir,paste0("occupancy_",Taxa,".rds"))
    
    geos <- length(unique(data$geo2))
    
    datMod <- data %>%
      tidyr::pivot_wider(names_from = "quart"
                         , values_from = "success"
                         ) %>%
      tidyr::nest(data = -c(year,geo2)) %>%
      dplyr::filter(geo2 == sort(geo2)[1]) %>%
      dplyr::mutate(umf = map(data
                              , function(x) unmarkedFrameOccu(y = x %>% dplyr::select(grep("q\\d{1}",names(x),value = TRUE)))
                              )
                  , mod = map(umf
                              , function(x) stan_occu(~1 ~1
                                                      , data = x
                                                      , chains = if(testing) testChains else useChains
                                                      , iter = if(testing) testIter else useIter
                                                      )
                              )
                  , occ = map(mod
                              , function(x) tibble(run = 1:draws
                                                   , occ = sample(plogis(extract(x,"beta_state[(Intercept)]")[[1]]),draws)
                                                   )
                              )
                  )
    
    
    
   
    
    write_rds(mod,outFile)
    
  }
  
  