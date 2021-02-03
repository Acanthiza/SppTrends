
  #---------Filter-------
  
  # Does an area in a year (for a taxonomic group) have enough lists?
  enoughLists <- dat %>%
    dplyr::distinct(Class,year,geo1,geo2,cell) %>%
    dplyr::count(Class,year,geo1,geo2,name="lists") %>%
    dplyr::filter(lists > 1)
  
  # Filter to areas*years*taxa that have enough lists
  datTemp <- dat %>%
    dplyr::inner_join(enoughLists) %>%
    dplyr::distinct(Class,Taxa,year,geo1,geo2)
  
  # Does a Taxa within an area have enough years where it appears on a list?
  enoughYears <- datTemp %>%
    dplyr::mutate(distinctYears = length(unique(.$year))
                  , thresh = 0.1*distinctYears
                  ) %>%
    dplyr::add_count(Class,Taxa,geo1,geo2, name = "years") %>%
    dplyr::filter(years > 0.2*distinctYears) %>%
    dplyr::distinct(Class,Taxa,geo1,geo2)
  
  # Does an area in a year have enough Taxa
  enoughTaxa <- datTemp %>%
    dplyr::inner_join(enoughYears) %>%
    dplyr::count(Class,geo1,geo2,year, name = "taxas") %>%
    dplyr::filter(taxas > 5) %>%
    dplyr::distinct(Class,geo1,geo2,year)
  
  datFilter <- dat %>%
    dplyr::inner_join(enoughYears) %>%
    dplyr::inner_join(enoughTaxa) %>%
    dplyr::filter(listLength > 1)
  
  # trials <- datFilter %>%
  #   dplyr::count(Class,year,geo,cell,name="listLength") %>%
  #   dplyr::count(Class,year,geo,name="trials")
  # 
  # success <- datFilter %>%
  #   dplyr::inner_join(trials) %>%
  #   dplyr::count(Class,Taxa,year,geo,name="success")
  # 
  taxaGeo <- datFilter %>%
    dplyr::distinct(Class,Taxa,geo1,geo2)
  
  
  #-------RR prep---------
  
  datRR <- datFilter %>%
    dplyr::group_by(Class,geo1,geo2,year) %>%
    dplyr::mutate(trials = n_distinct(cell)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Class,geo1,geo2,year,trials,Taxa) %>%
    dplyr::summarise(success = n()) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from = "Taxa", values_from = "success", values_fill = 0) %>%
    tidyr::pivot_longer((grep("trials",names(.))+1):ncol(.),names_to = "Taxa", values_to = "success") %>%
    dplyr::inner_join(taxaGeo) %>%
    dplyr::mutate(prop = success/trials)
  
  
  #------LL Prep---------------
  
  datLL <- datFilter %>%
    # dplyr::group_by(Class,geo1,geo2) %>%
    # dplyr::mutate(trials = n_distinct(cell)) %>%
    # dplyr::ungroup() %>%
    # dplyr::group_by(Class,geo1,geo2,year,trials,listLength,Taxa) %>%
    # dplyr::summarise(success = n()) %>%
    # dplyr::ungroup() %>%
    tidyr::nest(data = -c(Class,geo1,geo2,year
                          #,trials
                          ,listLength
                          )
                ) %>%
    dplyr::mutate(data = map(data
                             ,. %>%
                               tidyr::pivot_wider(names_from = "Taxa", values_from = "success", values_fill = 0) %>%
                               tidyr::pivot_longer(3:ncol(.),names_to = "Taxa", values_to = "success")
                             )
                  ) %>%
    tidyr::unnest(cols = c(data)) %>%
    dplyr::inner_join(taxaGeo)
  
  #--------RR explore------
  
  Y <- "prop"
  
  # variables to explore
  varExp <- c(get("Y")
              , colnames(datRR)
              ) %>%
    unique() %>%
    grep(pattern = "medYear",invert = TRUE, value=TRUE)
  
  datExp <- datRR %>%
    dplyr::select(all_of(varExp))
  
  #  Missing values
  plot_missing(datExp)
  
  # Count discrete variables
  ggplot(datExp %>%
           dplyr::mutate_if(is.factor,as.character) %>%
           dplyr::select_if(is.character) %>%
           tidyr::gather(variable,value,1:ncol(.)) %>%
           dplyr::group_by(variable) %>%
           dplyr::mutate(levels = n_distinct(value)) %>%
           dplyr::ungroup() %>%
           dplyr::filter(levels < maxLevels)
         ) +
    geom_histogram(aes(value),stat="count") +
    facet_wrap(~variable, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  # Count continuous variables
  ggplot(datExp %>%
           dplyr::select_if(is.numeric) %>%
           tidyr::gather(variable,value,1:ncol(.))
         , aes(value)
         ) +
    geom_histogram() +
    facet_wrap(~variable, scales = "free")
  
  # Y vs. Discrete
  ggplot(datExp %>%
           dplyr::mutate(UQ(rlang::sym(Y)) := factor(!!ensym(Y))) %>%
           dplyr::mutate_if(is.factor,as.character) %>%
           dplyr::select_if(is.character) %>%
           dplyr::mutate(UQ(rlang::sym(Y)) := as.numeric(!!ensym(Y))) %>%
           tidyr::gather(variable,value,2:ncol(.))
         ) +
    geom_boxplot(aes(value,!!ensym(Y))) +
    facet_wrap(~variable, scales = "free") +
    theme(axis.text.x=element_text(angle=90, vjust=0.5))
  
  # Y vs. Numeric
  ggplot(datExp %>%
           dplyr::select(varExp) %>%
           dplyr::select_if(is.numeric) %>%
           tidyr::gather(variable,value,2:ncol(.)) %>%
           dplyr::arrange(!!ensym(Y))
         , aes(value,!!ensym(Y))
         ) +
    geom_point(alpha = 0.5) +
    facet_wrap(~variable, scales = "free") +
    theme(axis.text.x=element_text(angle=90, vjust=0.5))
  
  #------RR model--------
  
  rr <- function(Taxa,rrDf) {
    
    print(Taxa)
    
    res <- list()
    
    geos <- length(unique(rrDf$geo2))
    
    if(geos > 1) {
      
      res$mod <- stan_glm(cbind(success,trials-success) ~ year*geo2
                      , data = rrDf
                      , family = binomial()
                      )
      
    } else {
      
      res$mod <- stan_glm(cbind(success,trials-success) ~ year
                        , data = rrDf
                        , family = binomial()
                        )
      
      }
    
    res$pred <- as_tibble(res$mod$data) %>%
      dplyr::distinct(year,geo1,geo2) %>%
      dplyr::mutate(col = row.names(.)
                    , success = 0
                    , trials = round(median(res$mod$data$trials),0)
                    ) %>%
      dplyr::left_join(as_tibble(posterior_predict(res$mod
                                                   , newdata = .
                                                   , re.form = NA#insight::find_formula(res$mod)$random[[1]]
                                                   )
                                 ) %>%
                         tibble::rownames_to_column(var = "row") %>%
                         tidyr::gather(col,value,2:ncol(.))
                       ) %>%
      dplyr::mutate(rawValue = as.numeric(value)
                    , value = rawValue/trials
                    )
    
    res$res <- res$pred %>%
      dplyr::group_by(year,geo1,geo2) %>%
      dplyr::summarise(modMean = mean(value)
                       , modMedian = quantile(value, 0.5)
                       , modci90lo = quantile(value, 0.05)
                       , modci90up = quantile(value, 0.95)
                       , text = paste0(round(modMedian,2)," (",round(modci90lo,2)," to ",round(modci90up,2),")")
                       ) %>%
      dplyr::ungroup()
    
    res$plot <- res$mod$data %>%
      dplyr::distinct(year,geo2) %>%
      tidybayes::add_fitted_draws(res$mod, n = 500) %>%
      ggplot(aes(x = year)) +
        geom_line(aes(y = .value, group = .draw), alpha = 0.05) +
        geom_line(data = res$res
                  , aes(y = modMean)
                  , size = 2
                  ) +
        geom_point(data = res$mod$data
                    , aes(y = prop
                          , colour = trials
                          )
                    ) +
        facet_grid(~geo2) +
        scale_colour_viridis_c() +
        labs(title = bquote(~italic(.(Taxa)))
             , subtitle = "Thick line is median credible value for that year"
             )
    
    write_rds(res,path(outDir,paste0("reporting-rate_",Taxa,".rds")))
    
    return(res)
    
  }
  
  
  ll <- function(Taxa,llDf) {
    
    print(Taxa)
    
    res <- list()
    
    geos <- length(unique(llDf$geo2))
    
    if(geos > 1) {
      
      res$mod <- stan_glm(success ~ year*geo2*log(listLength)
                          , data = llDf
                          , family = binomial()
                          )
      
    } else {
      
      res$mod <- stan_glm(success ~ year*log(listLength)
                          , data = llDf
                          , family = binomial()
                          )
      
    }
    
    res$pred <- res$mod$data %>%
      dplyr::distinct(geo1,geo2,year) %>%
      dplyr::mutate(listLength = median(llDf$listLength)
                    , col = row.names(.)
                    ) %>%
      dplyr::left_join(as_tibble(posterior_predict(res$mod
                                                   , newdata = .
                                                   , re.form = insight::find_formula(res$mod)$random
                                                   )
                                 ) %>%
                         tibble::rownames_to_column(var = "row") %>%
                         tidyr::gather(col,value,2:ncol(.))
                       )
    
    res$res <- res$pred %>%
      dplyr::group_by(year,geo1,geo2,listLength) %>%
      dplyr::summarise(n = n()
                       , nCheck = nrow(as_tibble(modLL))
                       , modMean = mean(value)
                       , modMedian = quantile(value, 0.5)
                       , modci90lo = quantile(value, 0.05)
                       , modci90up = quantile(value, 0.95)
                       , text = paste0(round(modMedian,2)," (",round(modci90lo,2)," to ",round(modci90up,2),")")
                       ) %>%
      dplyr::ungroup()
    
    res$plot <- res$mod$data %>%
      dplyr::distinct(year,geo2) %>%
      dplyr::full_join(tibble(probs = c(0.05
                                        , 0.5
                                        , 0.95
                                        )
                              ) %>%
                         dplyr::mutate(listLength = map_dbl(probs,~quantile(res$mod$data$listLength,probs = .)))
                       , by = character()
                       ) %>%
      dplyr::mutate(length = paste0("list length quantile ",probs," = ",listLength)) %>%
      tidybayes::add_fitted_draws(res$mod, n = 500) %>%
      ggplot(aes(x = year, y = success)) +
        geom_line(aes(y = .value, group = .draw), alpha = 0.05) +
        geom_jitter(data = res$mod$data
                    , aes(colour = listLength)
                    , height = 0.025
                    , width = 0.2
                    ) +
        facet_grid(length~geo2) +
        scale_colour_viridis_c() +
        labs(title = bquote(~italic(.(Taxa))))
    
    write_rds(res,path(outDir,paste0("list-length_",Taxa,".rds")))
    
    return(res)
    
  }
  
  
  taxaRes <- datRR %>%
    tidyr::nest(dataRR = c(geo1,geo2,year,success,trials,prop)) %>%
    dplyr::left_join(datLL %>%
                       tidyr::nest(dataLL = c(geo1,geo2,year,success,listLength,cell,list))
                     ) %>%
    dplyr::sample_n(5) %>% # TESTING
    dplyr::mutate(rr = future_map2(Taxa,dataRR,rr)
                  , ll = future_map2(Taxa,dataLL,ll)
                  )
  
  
  
  
  #------Residuals--------
  
  modResid <- tibble(residual = residuals(modRR)
                     , fitted = fitted(modRR)
                     ) %>%
    dplyr::bind_cols(as_tibble(modRR$data))
  
  
  # Y vs. residuals
  ggplot(modResid, aes(fitted,residual)) +
    geom_point(size = 2) +
    geom_hline(aes(yintercept = 0), linetype = 2, colour = "red") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) %>%
    scale_colour_viridis_d(end=0.9)
  
  # Time vs. residuals
  ggplot(modResid %>%
           dplyr::select_if(is.numeric) %>%
           tidyr::gather(variable,value,2:ncol(.))
         , aes(value,residual)
         ) +
    geom_point(size = 2) +
    geom_smooth(method = "lm") +
    facet_wrap(~variable
               , scales = "free_x"
               ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_colour_viridis_d()
  
  
  
    
  
  #-----TESTING--------
  
  test <- datLL %>%
    dplyr::filter(grepl(paste0(testSpp,collapse="|"),Taxa)) %>%
    dplyr::mutate(geo2 = fct_reorder(geo2,-success,mean))
  
  ggplot(test,aes(year
                  ,success
                  , colour = listLength
                  )
         ) +
    geom_jitter(height = 0.0005
                ,width = 0.005
                ,alpha = 0.5
                #,aes(colour=geo2)
                ) +
    geom_smooth() +
    facet_grid(~geo2)
  
  
  #-------LL model---------
  
  modLL <- stan_glm(success ~ year*geo2*log(listLength)
                    , data = test
                    , family = binomial()
                    )
  
  modPredLL <- as_tibble(modLL$data) %>%
    dplyr::distinct(Taxa,year,geo1,geo2) %>%
    dplyr::mutate(listLength = median(datLL$listLength)
                  , col = row.names(.)
                  ) %>%
    dplyr::left_join(as_tibble(posterior_predict(modLL
                                                 , newdata = .
                                                 , re.form = insight::find_formula(modLL)$random
                                                 )
                               ) %>%
                       tibble::rownames_to_column(var = "row") %>%
                       tidyr::gather(col,value,2:ncol(.))
                     )
  
  modResLL <- modPredLL %>%
    dplyr::group_by(Taxa,year,geo1,geo2,listLength) %>%
    dplyr::summarise(n = n()
                     , nCheck = nrow(as_tibble(modLL))
                     , modMean = mean(value)
                     , modMedian = quantile(value, 0.5)
                     , modci90lo = quantile(value, 0.05)
                     , modci90up = quantile(value, 0.95)
                     , text = paste0(round(modMedian,2)," (",round(modci90lo,2)," to ",round(modci90up,2),")")
                     ) %>%
    dplyr::ungroup()
  
  
  
  modLL$data %>%
    dplyr::distinct(Taxa,year,geo2) %>%
    dplyr::full_join(tibble(probs = c(0.05
                                      , 0.5
                                      , 0.95
                                      )
                            ) %>%
                       dplyr::mutate(listLength = map_dbl(probs,~quantile(modLL$data$listLength,probs = .)))
                     , by = character()
                     ) %>%
    dplyr::mutate(length = paste0("list length quantile ",probs," = ",listLength)) %>%
    tidybayes::add_fitted_draws(modLL, n = 200) %>%
    ggplot(aes(x = year, y = success)) +
      # tidybayes::stat_lineribbon(aes(y = .value)
      #                            , .width = c(.999, .95, .8, .5)
      #                            , alpha = 0.5
      #                            ) +
      geom_line(aes(y = .value, group = .draw), alpha = 0.05) +
      geom_jitter(data = modLL$data
                  , aes(colour = listLength)
                  , alpha = 0.01
                  , height = 0.025
                  , width = 0.2
                  ) +
      facet_grid(length~geo2)

    

  