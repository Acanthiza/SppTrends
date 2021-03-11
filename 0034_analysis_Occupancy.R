

  timer$start("occ")
  
  #-------Filter for analysis---------
  
  # Sampling unit for occupancy is based on a cell within a yearmon. Thus, success is presence within a cell in a yearmon
  
  allScales <- c("geo1","geo2","cell","site") #in order
  analysisScales <- allScales[-length(allScales)]
  
  success <- datTidy %>%
    dplyr::distinct(Taxa,year,yday,!!ensym(taxGroup),across(any_of(allScales))) %>%
    dplyr::mutate(success = 1)
  
  datFiltered <- success %>%
    dplyr::mutate(list = paste0(year,"-",yday,"-",!!ensym(taxGroup),"-",get(analysisScales[length(analysisScales)]))) %>%
    dplyr::add_count(list, name = "listLength") %>%
    filter_taxa_data(minListLengthThresh = 0
                     , maxListLengthOccurenceThresh = 0
                     , minlistOccurenceThresh = 0
                     , minYearsThresh = 3
                     , minListLengthsThresh = 0
                     , minCellsThresh = 3
                     , minSitesThresh = 3
                     , minYearSpanThresh = 10
                     #, allQuartVisits = TRUE
                     , timeVars = c("year","yday")
                     )
  
  
  #-------Data prep---------
  
  # data_pivot_fill_zero <- function(df) {
  #   
  #   df %>%
  #     tidyr::pivot_wider(names_from = c("Taxa","yday"), values_from = "success",values_fill = 0) %>%
  #     tidyr::pivot_longer(cols = contains(unique(taxaGeo$Taxa))
  #                         , names_to = c("Taxa","yday")
  #                         , names_sep = "_"
  #                         , values_to = "success"
  #                         )
  #   
  # }
  
  siteYearVisit <- datFiltered %>%
    dplyr::distinct(across(any_of(taxGroup)),geo1,geo2,cell,site,year,yday,list,listLength)
  
  taxaCellYear <- datFiltered %>%
    dplyr::distinct(across(any_of(taxGroup)),Taxa,cell,year)
  
  dat <- siteYearVisit %>%
    dplyr::left_join(taxaCellYear) %>%
    dplyr::left_join(datFiltered) %>%
    dplyr::mutate(success = if_else(is.na(success),0,success)) %>%
    purrr::when(testing ~ (.) %>% dplyr::filter(Taxa %in% tests)
                , !testing ~ (.)
                ) %>%
    tidyr::nest(data = -c(Order,Taxa)) %>%
    dplyr::left_join(luTax %>%
                       dplyr::select(Taxa,Common)
                     )
  
  #------model function--------
  
  occ <- function(Taxa,data,draws = 200, useGAM = FALSE) {
    
    print(paste0(Taxa))
    
    outFile <- fs::path(outDir,paste0("occupancy_",Taxa,".rds"))
    
    geos <- length(unique(data$geo2))
    
    datOcc <- data %>%
      dplyr::distinct(year,geo2,cell,yday,site,success) %>%
      tidyr::nest(data = -c(year,geo2,cell)) %>%
      dplyr::mutate(data = map(data, . %>%
                                 tidyr::pivot_wider(names_from = "yday"
                                                    , values_from = "success"
                                                    #, values_fill = 0
                                                    )
                               )
                    , trials = map_dbl(data,nrow)
                    ) %>%
      dplyr::filter(trials > 2) %>%
      
      #dplyr::sample_n(2) %>% # TESTING
      
      dplyr::mutate(umf = map(data
                              , function(x) unmarkedFrameOccu(y = x %>% dplyr::select(-1))
                              )
                  , modYear = map(umf
                              , function(x) unmarked::occu(~1 ~1
                                                           , data = x
                                                           )
                              )
                  , occ = map_dbl(modYear
                              , ~backTransform(.,type = "state")@estimate
                              )
                   , occ = if_else(occ == 0, 0.000000000001,occ)
                   , occ = if_else(occ == 1, 0.999999999999,occ)
                  , det = map_dbl(modYear
                              , ~backTransform(.,type = "det")@estimate
                              )
                  ) %>%
      dplyr::select(where(Negate(is.list))) %>%
      dplyr::mutate(across(where(is.character),factor))
    
    write_feather(datOcc,path(outDir,paste0("occupancyDf_",Taxa,".feather")))
    
    # GAM
    if(useGAM) {
      
      if(geos > 1) {
        
        mod <- stan_gamm4(occ ~ s(year, k = 4, bs = "ts") + geo2 + s(year, k = 4, by = geo2, bs = "ts")
                          , data = datOcc
                          , family = binomial
                          , chains = if(testing) testChains else useChains
                          , iter = if(testing) testIter else useIter
                          )
        
      } else {
        
        mod <- stan_gamm4(occ ~ s(year, k = 4)
                          , data = datOcc
                          , family = binomial
                          , chains = if(testing) testChains else useChains
                          , iter = if(testing) testIter else useIter
                          )
        
      }
      
    } else {
      
      # GLM
      if(geos > 1) {

        mod <- stan_betareg(occ ~ year*geo2
                            , data = datOcc
                            , chains = if(testing) testChains else useChains
                            , iter = if(testing) testIter else useIter
                            )

      } else {

        mod <- stan_betareg(occ ~ year
                            , data = datOcc
                            , chains = if(testing) testChains else useChains
                            , iter = if(testing) testIter else useIter
                            )

      }
      
      
    }
    
    write_rds(mod,outFile)
    
  }
  
  #------Run models-------
  
  todo <- dat %>%
    dplyr::mutate(outFile = fs::path(outDir,paste0("occupancy_",Taxa,".rds"))
                  , done = map_lgl(outFile,file.exists)
                  ) %>%
    dplyr::filter(!done)
  
  if(nrow(todo) > 0) {
    
    if(nrow(todo) > useCores/(if(testing) testChains else useChains)) {
      
      future_pwalk(list(todo$Taxa
                        , todo$data
                        )
                   , occ
                   , useGAM = FALSE
                   )
      
    } else {
      
      pwalk(list(todo$Taxa,todo$data),occ)
      
    }
    
    
  }
  
  
  #--------Explore models-----------
  
  taxaModsOC <- dat %>%
    dplyr::mutate(mod = fs::path(outDir,paste0("occupancy_",Taxa,".rds"))
                  , data = fs::path(outDir,paste0("occupancyDf_",Taxa,".feather"))
                  , exists = map_lgl(mod,file.exists)
                  ) %>%
    dplyr::filter(exists) %>%
    dplyr::mutate(data = map(data,read_feather)
                  , mod = map(mod,read_rds)
                  , type = "Occupancy"
                  , res = pmap(list(Taxa
                                  , Common
                                  , data
                                  , mod
                                  , type
                                  , respVar = "occ"
                                  )
                             , mod_explore
                             )
                  )
  
  timer$stop("occ", comment = paste0("Occupancy models run for "
                                     ,nrow(taxaMods)
                                     ," taxa, of which "
                                     ,nrow(todo)
                                     ," were new"
                                     )
             )
  
  
  