
  timer$start("rr")
  
  #-------Filter for analysis---------
  
  rrGroup <- alwaysGroup %>% grep("yday|month",.,invert = TRUE, value = TRUE)
  
  datFiltered <- datTidy %>%
    dplyr::distinct(Taxa,across(all_of(taxGroup)),across(all_of(rrGroup))) %>%
    add_list_length(groups = c(taxGroup,rrGroup)) %>%
    filter_taxa_data()
  
  
  #------Absences from cooccurs-------
  
  datCooccur <- datFiltered %>%
    dplyr::mutate(p = 1) %>%
    dplyr::bind_rows(datFiltered %>%
                       dplyr::distinct(!!ensym(taxGroup)
                                      , geo2
                                      ) %>%
                       dplyr::inner_join(contextCooccur) %>%
                       tidyr::unnest(cols = c(pair)) %>%
                       dplyr::select(where(negate(is.list))) %>%
                       dplyr::rename(indicatesAbsence = Taxa
                                     , Taxa = sp2
                                     ) %>%
                       dplyr::inner_join(datFiltered) %>%
                       dplyr::mutate(p = 0) %>%
                       dplyr::select(-Taxa) %>%
                       dplyr::rename(Taxa = indicatesAbsence)
                     ) %>%
    dplyr::group_by(across(all_of(grep("^p$",names(datFiltered),value = TRUE,invert = TRUE)))) %>%
    dplyr::summarise(p = max(p)) %>%
    dplyr::ungroup()
  
  
  #-------Data prep---------
  
  dat <- datCooccur %>%
    dplyr::group_by(Taxa,!!ensym(taxGroup),across(all_of(rrGroup[-length(rrGroup)]))) %>%
    dplyr::summarise(success = sum(p)
                     , trials = n()
                     ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(prop = success/trials) %>%
    tidyr::nest(data = -c(!!ensym(taxGroup),Taxa)) %>%
    dplyr::left_join(luTax %>%
                       dplyr::select(Taxa,Common)
                     )
  
  
  #------Model function--------
  
  rr <- function(Taxa,data) {
    
    print(Taxa)
    
    outFile <- fs::path(outDir,paste0("reporting-rate_Mod_",Taxa,".rds"))
    
    geos <- length(unique(data$geo2))
    
    cells <- length(unique(data$cell))
    
    # GAM
    if(geos > 1) {
      
      mod <- stan_gamm4(cbind(success,trials - success) ~
                          s(year, k = 4, bs = "ts") +
                          geo2 +
                          s(year, k = 4, by = geo2, bs = "ts")
                        , data = data
                        , family = binomial()
                        , random = ~(1|cell)
                        , chains = if(testing) testChains else useChains
                        , iter = if(testing) testIter else useIter
                        )
      
    } else {
      
      mod <- stan_gamm4(cbind(success,trials-success) ~ s(year, k = 4, bs = "ts")
                            , data = data
                            , random = if(cells > 1) formula(~ (1|cell)) else NULL
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
    {if(testing) (.) %>% dplyr::filter(Taxa %in% tests) else (.)} %>%
    dplyr::mutate(outFile = fs::path(outDir,paste0("reporting-rate_Mod_",Taxa,".rds"))
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
  
  taxaModsRR <- dat %>%
    dplyr::mutate(modPath = fs::path(outDir,paste0("reporting-rate_Mod_",Taxa,".rds"))) %>%
    dplyr::filter(file.exists(modPath)) %>%
    dplyr::mutate(type = "Reporting rate")
  
  doof <- taxaModsRR
  
  future_pwalk(list(doof$Taxa
                    , doof$Common
                    , doof$data
                    , doof$modPath
                    , doof$type
                    )
               , mod_explore
               , respVar = "prop"
               )
   
  timer$stop("rr", comment = paste0("Reporting rate models run for "
                                    ,nrow(taxaModsRR)
                                    ," taxa, of which "
                                    ,nrow(todo)
                                    ," were new"
                                    )
             )
  