
  timer$start("rr")
  
  #-------Filter for analysis---------
  
  # Sampling unit, or list, for rr and ll is based on a cell within a year. Thus, success is presence within a cell in a year
  
  allScales <- c("geo1","geo2","cell","site") #in order
  analysisScales <- allScales[-length(allScales)]
  
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
    filter_taxa_data()
  
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
    dplyr::mutate(prop = success/trials) %>%
    tidyr::nest(data = c(any_of(analysisScales),year,list,listLength,success,trials,prop)) %>%
    dplyr::left_join(luTax %>%
                       dplyr::select(Taxa,Common)
                     )
  
  
  #------Model function--------
  
  rr <- function(Taxa,data) {
    
    print(Taxa)
    
    outFile <- fs::path(outDir,paste0("reporting-rate_",Taxa,".rds"))
    
    geos <- length(unique(data$geo2))
    
    # GAM
    if(geos > 1) {
      
      mod <- stan_gamm4(cbind(success,trials - success) ~ s(year, k = 4) + geo2 + s(year, k = 4, by = geo2)
                            , data = data
                            , family = binomial()
                            , random = ~(1|cell)
                            , chains = if(testing) testChains else useChains
                            , iter = if(testing) testIter else useIter
                            )
      
    } else {
      
      mod <- stan_gamm4(cbind(success,trials-success) ~ s(year, k = 4)
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
 
  
  # Check if rr models have been run - run if not
  todo <- dat %>%
    dplyr::mutate(outFile = fs::path(outDir,paste0("reporting-rate_",Taxa,".rds"))
                  , done = map_lgl(outFile,file.exists)
                  ) %>%
    dplyr::filter(!done)
  
  if(nrow(todo) > 0) {
    
    if(nrow(todo) > useCores/(if(testing) testChains else useChains)) {
      
      future_pwalk(list(todo$Taxa
                        , todo$data
                        )
                   ,rr
                   )
      
    } else {
      
      pwalk(list(todo$Taxa,todo$data),rr)
      
    }
    
    
  }
  
  
  #--------Explore models-----------
  
  taxaMods <- dat %>%
    dplyr::mutate(mod = fs::path(outDir,paste0("reporting-rate_",Taxa,".rds"))
                  , exists = map_lgl(mod,file.exists)
                  ) %>%
    dplyr::filter(exists) %>%
    dplyr::mutate(mod = map(mod,read_rds)) %>%
    dplyr::mutate(res = pmap(list(Taxa
                                  , Common
                                  , data
                                  , mod
                                  , modType = "Reporting rate"
                                  )
                             , mod_explore
                             )
                  )
   
  timer$stop("rr", comment = paste0("Reporting rate models run for "
                                    ,nrow(taxaMods)
                                    ," taxa, of which "
                                    ,nrow(todo)
                                    ," were new"
                                    )
             )
  