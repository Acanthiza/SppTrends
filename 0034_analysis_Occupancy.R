

  timer$start("occ")
  
  #-------Filter for analysis---------
  
  # Sampling unit for occupancy is based on a cell within a yearmon. Thus, success is presence within a cell in a yearmon
  
  occGroup <- alwaysGroup
  
  datFiltered <- datTidy %>%
    dplyr::distinct(Taxa,across(all_of(taxGroup)),across(all_of(occGroup))) %>%
    add_list_length(groups = c(taxGroup,occGroup)) %>%
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
    tidyr::nest(data = -c(!!ensym(taxGroup),Taxa)) %>%
    dplyr::left_join(luTax %>%
                       dplyr::select(Taxa,Common)
                     )
  
  #------model function--------
  
  occ <- function(Taxa,data,draws = 200, useGAM = TRUE) {
    
    print(paste0(Taxa))
    
    outFile <- fs::path(outDir,paste0("occupancy_Mod_",Taxa,".rds"))
    
    geos <- length(unique(data$geo2))
    
    cells <- length(unique(data$cell))
    
    datOccPrep <- data %>%
      dplyr::distinct(year,geo2,cell,yday,site,p) %>%
      tidyr::nest(data = -c(year,geo2,cell)) %>%
      dplyr::filter(map_lgl(data,~sum(.$p) > 0)) %>%
      dplyr::mutate(data = map(data, . %>%
                                 tidyr::pivot_wider(names_from = "yday"
                                                    , values_from = "p"
                                                    #, values_fill = 0
                                                    )
                               )
                    , trials = map_dbl(data,nrow)
                    ) %>%
      dplyr::filter(trials > 2) %>%
      dplyr::group_by(geo2) %>%
      dplyr::mutate(years = n_distinct(year)) %>%
      dplyr::filter(years > 3)
    
    if(nrow(datOccPrep) > 3) {
      
      datOcc <- datOccPrep %>%
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
                      , occ = if_else(occ == 0, 0.00000001,occ)
                      , occ = if_else(occ == 1, 0.99999999,occ)
                    , det = map_dbl(modYear
                                , ~backTransform(.,type = "det")@estimate
                                )
                    ) %>%
        dplyr::select(where(Negate(is.list))) %>%
        dplyr::group_by(geo2) %>%
        dplyr::mutate(years = n_distinct(year)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(years > 3) %>%
        dplyr::mutate(across(where(is.character),factor))
      
        write_feather(datOcc,path(outDir,paste0("occupancyDf_",Taxa,".feather")))
        
      geos <- length(unique(datOcc$geo2))
      randCell <- length(unique(datOcc$cell)) > 1 & min(table(datOcc$cell)) > 1
      
        # GAM
      if(useGAM) {
        
        if(geos > 1) {
          
          mod <- stan_gamm4(occ ~ s(year, k = 4, bs = "ts") + geo2 + s(year, k = 4, by = geo2, bs = "ts")
                            , data = datOcc
                            , family = mgcv::betar()
                            , random = if(randCell) formula(~ (1|cell)) else NULL
                            , chains = if(testing) testChains else useChains
                            , iter = if(testing) testIter else useIter
                            )
          
        } else {
          
          mod <- stan_gamm4(occ ~ s(year, k = 4)
                            , data = datOcc
                            , family = mgcv::betar()
                            , random = if(randCell) formula(~ (1|cell)) else NULL
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
    
  }
  
  #------Run models-------
  
  todo <- dat %>%
    {if(testing) (.) %>% dplyr::filter(Taxa %in% tests) else (.)} %>%
    dplyr::mutate(outFile = fs::path(outDir,paste0("occupancy_Mod_",Taxa,".rds"))
                  , done = map_lgl(outFile,file.exists)
                  ) %>%
    dplyr::filter(!done)
  
  if(nrow(todo) > 0) {
    
    if(nrow(todo) > useCores/(if(testing) testChains else useChains)) {
      
      future_pwalk(list(todo$Taxa
                        , todo$data
                        )
                   , occ
                   )
      
    } else {
      
      pwalk(list(todo$Taxa,todo$data),occ)
      
    }
    
    
  }
  
  
  #--------Explore models-----------
  
 taxaModsOC <- dat %>%
    dplyr::mutate(data = fs::path(outDir,paste0("occupancyDf_",Taxa,".feather"))
                  , modPath = fs::path(outDir,paste0("occupancy_Mod_",Taxa,".rds"))
                  ) %>%
    dplyr::filter(file.exists(modPath)) %>%
    dplyr::mutate(data = map(data,read_feather)
                  , type = "Occupancy"
                  )
  
  doof <- taxaModsOC
  
  future_pwalk(list(doof$Taxa
                    , doof$Common
                    , doof$data
                    , doof$modPath
                    , doof$type
                    )
               , mod_explore
               , respVar = "occ"
               )
  
  timer$stop("occ", comment = paste0("Occupancy models run for "
                                     ,nrow(taxaModsOC)
                                     ," taxa, of which "
                                     ,nrow(todo)
                                     ," were new"
                                     )
             )
  
  
  