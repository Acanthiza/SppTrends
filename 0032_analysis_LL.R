
  timer$start("ll")
  
  #-------Filter for analysis---------
  
  # Sampling unit, or list, for rr and ll is based on a cell within a year. Thus, success is presence within a cell in a year
  
  allScales <- c("geo1","geo2","cell") #in order
  analysisScales <- allScales
  
  trials <- datTidy %>%
    dplyr::distinct(year,!!ensym(taxGroup),across(any_of(allScales))) %>%
    dplyr::group_by(year,!!ensym(taxGroup),across(any_of(analysisScales))) %>%
    dplyr::summarise(trials = n()) %>%
    dplyr::ungroup()
  
  success <- datTidy %>%
    dplyr::distinct(Taxa,year,!!ensym(taxGroup),across(any_of(allScales))) %>%
    dplyr::group_by(Taxa,year,!!ensym(taxGroup),across(any_of(analysisScales))) %>%
    dplyr::summarise(success = n()) %>%
    dplyr::ungroup()
  
  datFiltered <- trials %>%
    dplyr::left_join(success) %>%
    dplyr::mutate(list = paste0(year,"-",!!ensym(taxGroup),"-",get(analysisScales[length(analysisScales)]))) %>%
    dplyr::add_count(list, name = "listLength") %>%
    filter_taxa_data(minCellsThresh = 0)
  
  # trials <- datTidy %>%
  #   dplyr::inner_join(datFiltered) %>%
  #   dplyr::distinct(year,!!ensym(taxGroup),across(any_of(allScales))) %>%
  #   dplyr::group_by(year,!!ensym(taxGroup),across(any_of(analysisScales))) %>%
  #   dplyr::summarise(trials = n()) %>%
  #   dplyr::ungroup()
  # 
  # success <- datFiltered %>%
  #   dplyr::distinct(Taxa,year,!!ensym(taxGroup),across(any_of(allScales))) %>%
  #   dplyr::group_by(Taxa,year,!!ensym(taxGroup),across(any_of(analysisScales))) %>%
  #   dplyr::summarise(success = n()) %>%
  #   dplyr::ungroup()
  # 
  # datForAnalysis <- trials %>%
  #   dplyr::left_join(success) %>%
  #   dplyr::mutate(list = paste0(year,"-",!!ensym(taxGroup),"-",get(analysisScales[length(analysisScales)]))) %>%
  #   dplyr::add_count(list, name = "listLength")
  
  taxaGeo <- datFiltered %>%
    dplyr::distinct(!!ensym(taxGroup)
                    ,Taxa
                    ,across(any_of(analysisScales))
                    )
  
  #-------Data prep---------
  
  dat <- datFiltered %>%
    dplyr::inner_join(taxaGeo) %>%
    tidyr::pivot_wider(names_from = "Taxa", values_from = "success", values_fill = 0) %>%
    tidyr::pivot_longer(cols = any_of(unique(taxaGeo$Taxa)),names_to = "Taxa", values_to = "success") %>%
    dplyr::inner_join(taxaGeo) %>%
    purrr::when(testing ~ (.) %>% dplyr::filter(Taxa %in% tests)
                , !testing ~ (.)
                ) %>%
    dplyr::mutate(prop = success/trials
                  , listLengthLog = log(listLength)
                  ) %>%
    tidyr::nest(data = c(any_of(analysisScales),year,list,contains("listLength"),success,trials,prop)) %>%
    dplyr::left_join(luTax %>%
                       dplyr::select(Taxa,Common)
                     )
  
  
  
  #------model function--------
  
  ll <- function(Taxa,data) {
    
    print(Taxa)
    
    outFile <- fs::path(outDir,paste0("list-length_",Taxa,".rds"))
    
    geos <- length(unique(data$geo2))
    
    # GAM
    if(geos > 1) {
      
      mod <- stan_gamm4(cbind(success,trials - success) ~
                          s(year, k = 4, bs = "ts") +
                          s(year, k = 4, by = geo2, bs = "ts") +
                          s(year, k = 4, by = listLengthLog, bs = "ts") +
                          geo2 +
                          geo2*listLengthLog
                        , data = data
                        , family = binomial()
                        , random = ~(1|cell)
                        , chains = if(testing) testChains else useChains
                        , iter = if(testing) testIter else useIter
                        )
      
    } else {
      
      mod <- stan_gamm4(cbind(success,trials-success) ~ s(year, k = 4, bs = "ts") +
                          s(year, k = 4, by = listLengthLog, bs = "ts") +
                          listLengthLog
                        , data = data
                        , random = ~ (1|cell)
                        , family = binomial()
                        , chains = if(testing) testChains else useChains
                        , iter = if(testing) testIter else useIter
                        )
      
    }
    
    write_rds(mod,outFile)
    
  }
  
  
  #------Run models-------
 
  todo <- dat %>%
    dplyr::mutate(outFile = fs::path(outDir,paste0("list-length_",Taxa,".rds"))
                  , done = map_lgl(outFile,file.exists)
                  ) %>%
    dplyr::filter(!done)
  
  if(nrow(todo) > 0) {
    
    if(nrow(todo) > useCores/(if(testing) testChains else useChains)) {
      
      future_pwalk(list(todo$Taxa
                      , todo$data
                      )
               ,ll
               )
      
    } else {
      
      pwalk(list(todo$Taxa,todo$data),ll)
      
    }
    
    
  }
  
  
  #--------Explore models-----------
  taxaModsLL <- dat %>%
    dplyr::mutate(mod = fs::path(outDir,paste0("list-length_",Taxa,".rds"))
                  , exists = map_lgl(mod,file.exists)
                  ) %>%
    dplyr::filter(exists) %>%
    dplyr::mutate(mod = map(mod,read_rds)
                  , type = "List length"
                  , res = pmap(list(Taxa
                                  , Common
                                  , data
                                  , mod
                                  , type
                                  )
                             , mod_explore
                             )
                  )
  
  
  timer$stop("ll", comment = paste0("Reporting rate models run for "
                                    ,nrow(taxaModsLL)
                                    ," taxa, of which "
                                    ,nrow(todo)
                                    ," were new"
                                    )
             )
  
  