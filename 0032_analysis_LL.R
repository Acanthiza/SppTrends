
  timer$start("ll")
  
  #-------Filter for analysis---------
  
  # Sampling unit, or list, for rr and ll is based on a cell within a year. Thus, success is presence within a cell in a year
  
  llGroup <- alwaysGroup %>% grep("yday|month|site",.,invert = TRUE, value = TRUE)
  
  datFiltered <- datTidy %>%
    dplyr::distinct(Taxa,across(all_of(taxGroup)),across(all_of(llGroup)),site) %>%
    add_list_length(groups = c(taxGroup,llGroup)) %>%
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
    dplyr::summarise(p = sum(p)) %>%
    dplyr::ungroup()
  
  
  #-------Data prep---------
  
  dat <- datCooccur %>%
    dplyr::group_by(across(all_of(grep("site|^p$",names(datCooccur),value = TRUE,invert = TRUE)))) %>%
    dplyr::summarise(success = sum(p)
                     , trials = n()
                     ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(prop = success/trials
                  , listLengthLog = log(listLength)
                  ) %>%
    tidyr::nest(data = -c(!!ensym(taxGroup),Taxa)) %>%
    dplyr::left_join(luTax %>%
                       dplyr::select(Taxa,Common)
                     )
  
  
  #------model function--------
  
  ll <- function(Taxa,data) {
    
    print(Taxa)
    
    outFile <- fs::path(outDir,paste0("list-length_Mod_",Taxa,".rds"))
    
    geos <- length(unique(data$geo2))
    
    randCell <- length(unique(data$cell)) > 1 & max(table(factor(data$cell))) > 1
    
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
                        , random = if(randCell) formula(~ (1|cell)) else NULL
                        , chains = if(testing) testChains else useChains
                        , iter = if(testing) testIter else useIter
                        )
      
    } else {
      
      mod <- stan_gamm4(cbind(success,trials-success) ~ s(year, k = 4, bs = "ts") +
                          s(year, k = 4, by = listLengthLog, bs = "ts") +
                          listLengthLog
                        , data = data
                        , random = if(randCell) formula(~ (1|cell)) else NULL
                        , family = binomial()
                        , chains = if(testing) testChains else useChains
                        , iter = if(testing) testIter else useIter
                        )
      
    }
    
    write_rds(mod,outFile)
    
  }
  
  
  #------Run models-------
 
  todo <- dat %>%
    {if(testing) (.) %>% dplyr::filter(Taxa %in% tests) else (.)} %>%
    dplyr::mutate(outFile = fs::path(outDir,paste0("list-length_Mod_",Taxa,".rds"))
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
    dplyr::mutate(modPath = fs::path(outDir,paste0("list-length_Mod_",Taxa,".rds"))) %>%
    dplyr::filter(file.exists(modPath)) %>%
    dplyr::mutate(type = "List length")
  
  doof <- taxaModsLL
  
  future_pwalk(list(doof$Taxa
                    , doof$Common
                    , doof$data
                    , doof$modPath
                    , doof$type
                    )
               , mod_explore
               , respVar = "prop"
               )
  
  
  timer$stop("ll", comment = paste0("List length models run for "
                                    ,nrow(taxaModsLL)
                                    ," taxa, of which "
                                    ,nrow(todo)
                                    ," were new"
                                    )
             )
  
  