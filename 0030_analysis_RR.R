
  timer$start("rr")
  
  #-------Filter for analysis---------
  
  rrGroup <- alwaysGroup %>% grep("yday|month|site",.,invert = TRUE, value = TRUE)
  
  datLong <- datTidy %>%
    dplyr::distinct(Taxa,across(rrGroup)) %>%
    dplyr::count(Taxa,geo2)
  
  datWide <- datLong %>%
    tidyr::pivot_wider(names_from = "geo2", values_from = "n", values_fill = 0)
  
  dfWide <- data.frame(datWide[,-1])
  
  rownames(dfWide) <- datWide$Taxa
  
  datDist <- vegan::vegdist(dfWide
                   , method = "bray"
                   , binary = TRUE
                   )
    
  datDend <- clustMethod %>%
    dplyr::mutate(dend = map(method
                             ,~fastcluster::hclust(datDist, .)
                             )
                  )
  
  make_sil <- function(clustDf, distObj = datDist, clustID = "clust"){
    
    clustCol <- if(is.character(clustID)) which(names(clustDf)==clustID) else siteID
    
    silhouette(dplyr::pull(clustDf,clustCol),distObj)
    
  }
  
  # Turn an object of class silhouette into a data frame with one row per site
  make_sil_df <- function(clustDf,silObj) {
    
    clustDf %>%
      dplyr::bind_cols(tibble(neighbour = silObj[,2],sil_width = silObj[,3]))
    
  }
  
  # UPTOHERE - this is working, but not sure it is what is needed. Also need to get 'recGroup' by IBRA Subregion (geo2)
  datClust <- datDend %>%
    dplyr::mutate(clusters = map(dend
                                 , cutree
                                 , possibleGroups[possibleGroups < 0.5*nrow(datWide)]
                                 )
                  , clusters = map(clusters,as_tibble)
                  ) %>%
    dplyr::select(-dend) %>%
    tidyr::unnest(clusters) %>%
    tidyr::gather(groups,clust,2:ncol(.)) %>%
    dplyr::mutate(groups = as.integer(groups)) %>%
    tidyr::nest(clusters = c(clust)) %>%
    dplyr::mutate(clusters = map(clusters, . %>%
                                   dplyr::bind_cols(datWide[,1])
                                 )
                  ) %>%
    #dplyr::sample_n(5) %>% # TESTING
    dplyr::mutate(sil = map(clusters,make_sil)
                  , silDf = map2(clusters,sil,make_sil_df)
                  , macroSil = map_dbl(silDf,~mean(.$sil_width))
                  ) %>%
    tidyr::unnest(cols = c(silDf)) %>%
    dplyr::group_by(Taxa) %>%
    dplyr::filter(sil_width == max(sil_width)) %>%
    dplyr::filter(macroSil == max(macroSil)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(recGroup = map2(clusters
                                  ,clust
                                  ,function(x,y) x %>%
                                    dplyr::filter(clust == y) %>%
                                    dplyr::select(recGroup = Taxa)
                                  )
                  ) %>%
    tidyr::unnest(cols = c(recGroup)) %>%
    dplyr::select(where(negate(is.list)))
  
  plot(dend)
  
  
    filter_taxa_data()
  
  taxaCell <- datFiltered %>%
    dplyr::distinct(Taxa,cell)
  
  #-------Data prep---------
  
  
  
  dat <- datFiltered %>%
    dplyr::mutate(success = 1) %>%
    tidyr::nest(data = -c(!!ensym(taxGroup))) %>%
    dplyr::mutate(data = map(data
                             , . %>%
                               tidyr::pivot_wider(names_from = "Taxa", values_from = "success", values_fill = 0) %>%
                               tidyr::pivot_longer(cols = any_of(unique(taxaGeo$Taxa)),names_to = "Taxa", values_to = "success")
                             )
                  ) %>%
    tidyr::unnest(cols = c(data)) %>%
    dplyr::inner_join(taxaCell) %>%
    purrr::when(testing ~ (.) %>% dplyr::filter(Taxa %in% tests)
                , !testing ~ (.)
                ) %>%
    dplyr::group_by(Taxa,across(taxGroup),across(alwaysGroup %>% grep("site",.,invert = TRUE, value = TRUE))) %>%
    dplyr::summarise(success = sum(success)
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
  