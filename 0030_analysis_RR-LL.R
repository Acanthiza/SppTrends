
  timer$start("rrll")
  
  #-------datForRR---------
  
  # Sampling unit, or list, for rr and ll is based on a cell within a year. Thus, success is presence within a cell in a year
  
  useScale <- c("geo1","geo2","cell") #in order
  
  trials <- dat %>%
    dplyr::distinct(year,!!ensym(taxGroup),across(any_of(useScale))) %>%
    dplyr::group_by(year,!!ensym(taxGroup),across(any_of(useScale))) %>%
    dplyr::summarise(trials = n()) %>%
    dplyr::ungroup()
  
  success <- dat %>%
    dplyr::distinct(Taxa,year,!!ensym(taxGroup),across(any_of(useScale))) %>%
    dplyr::group_by(Taxa,year,!!ensym(taxGroup),across(any_of(useScale))) %>%
    dplyr::summarise(success = n()) %>%
    dplyr::ungroup()
  
  datFiltered <- trials %>%
    dplyr::left_join(success) %>%
    dplyr::mutate(list = paste0(year,"-",!!ensym(taxGroup),"-",useScale[length(useScale)])) %>%
    dplyr::add_count(list, name = "listLength") %>%
    filter_taxa_data()
  
  trials <- datFiltered %>%
    dplyr::distinct(year,!!ensym(taxGroup),across(any_of(useScale))) %>%
    dplyr::group_by(year,!!ensym(taxGroup),across(any_of(useScale))) %>%
    dplyr::summarise(trials = n()) %>%
    dplyr::ungroup()
  
  success <- datFiltered %>%
    dplyr::distinct(Taxa,year,!!ensym(taxGroup),across(any_of(useScale))) %>%
    dplyr::group_by(Taxa,year,!!ensym(taxGroup),across(any_of(useScale))) %>%
    dplyr::summarise(success = n()) %>%
    dplyr::ungroup()
  
  datForRR <- trials %>%
    dplyr::left_join(success) %>%
    dplyr::mutate(list = paste0(year,"-",!!ensym(taxGroup),"-",UQ(rlang::sym(useScale[length(useScale)])))) %>%
    dplyr::add_count(list, name = "listLength")
  
  taxaGeo <- datForRR %>%
    dplyr::distinct(!!ensym(taxGroup)
                    ,Taxa
                    ,across(any_of(useScale))
                    ) 
  
  #-------RR prep---------
  
  datRR <- datForRR %>%
    dplyr::inner_join(taxaGeo) %>%
    tidyr::pivot_wider(names_from = "Taxa", values_from = "success", values_fill = 0) %>%
    tidyr::pivot_longer(cols = any_of(unique(taxaGeo$Taxa)),names_to = "Taxa", values_to = "success") %>%
    dplyr::inner_join(taxaGeo) %>%
    purrr::when(testing ~ (.) %>% dplyr::filter(Taxa %in% tests)
                , !testing ~ (.)
                ) %>%
    tidyr::nest(data = c(geo1,geo2,year,cell,list,listLength,success,trials)) %>%
    dplyr::left_join(luTax %>%
                       dplyr::select(Taxa,Common)
                     )
  
  
  #------Functions--------
  
  explore <- function(Taxa,Common,df
                      ,respVar = "prop"
                      ,expVar = c("geo1","geo2","year")
                      ,maxLevels = 30
                      ) {
    
    res <- list()
    
    plotTitles <- bquote(~italic(.(Taxa))*":" ~ .(Common))
    
    if(!respVar %in% names(df)) df <- df %>%
      dplyr::group_by(across(any_of(expVar))) %>%
      dplyr::summarise(!!ensym(respVar) := sum(success)/n()) %>%
      dplyr::ungroup()
    
    Y <- respVar
    
    # variables to explore
    varExp <- c(respVar
                , colnames(df)
                ) %>%
      unique()
    
    datExp <- df %>%
      dplyr::select(all_of(varExp))
    
    hasNumeric <- datExp %>%
      dplyr::select(-1) %>%
      dplyr::select(where(is.numeric)) %>%
      ncol() %>%
      `>` (0)
    
    hasCharacter <- datExp %>%
      dplyr::select(-1) %>%
      dplyr::mutate(across(where(is.factor),as.character)) %>%
      dplyr::select(where(is.character)) %>%
      ncol() %>%
      `>` (0)
    
    #  Missing values
    res$missing <- plot_missing(datExp)
    
    # Character variables
    if(hasCharacter) {
      
      # count character
      res$countChar <- ggplot(datExp %>%
                              dplyr::mutate(across(where(is.factor),as.character)) %>%
                              dplyr::select_if(is.character) %>%
                              tidyr::gather(variable,value,1:ncol(.)) %>%
                              dplyr::group_by(variable) %>%
                              dplyr::mutate(levels = n_distinct(value)) %>%
                              dplyr::ungroup() %>%
                              dplyr::filter(levels < maxLevels)
                            ) +
        geom_histogram(aes(value)
                       , stat = "count"
                       ) +
        facet_wrap(~variable, scales = "free") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        labs(title = plotTitles
             , subtitle = "Count of levels within character variables"
             )
      
      # Y vs character
      res$YvsChar <-ggplot(datExp %>%
                             dplyr::mutate({{Y}} := factor(!!ensym(Y))) %>%
                             dplyr::mutate_if(is.factor,as.character) %>%
                             dplyr::select_if(is.character) %>%
                             dplyr::mutate({{Y}} := as.numeric(!!ensym(Y))) %>%
                             tidyr::gather(variable,value,2:ncol(.))
                           ) +
        geom_boxplot(aes(value,!!ensym(Y))) +
        facet_wrap(~variable, scales = "free") +
        theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
        labs(title = plotTitles
             , subtitle = paste0("Boxplots of response variable (",Y,") against character variables")
             )
      
    }
    
    # Numeric variables
    if(hasNumeric) {
      
      # Count numeric
      res$countNum <- ggplot(datExp %>%
                               dplyr::select(where(is.numeric)) %>%
                               tidyr::gather(variable,value,1:ncol(.))
                             , aes(value)
                             ) +
        geom_histogram() +
        facet_wrap(~variable, scales = "free") +
        labs(title = plotTitles
             , subtitle = "Histograms of numeric variables"
             )
      
      # Y vs. Numeric
      res$YvsNum <- ggplot(datExp %>%
                             dplyr::select(any_of(varExp)) %>%
                             dplyr::select(where(is.numeric)) %>%
                             tidyr::gather(variable,value,2:ncol(.)) %>%
                             dplyr::arrange(!!ensym(Y))
                           , aes(value,!!ensym(Y))
                           ) +
        geom_point(alpha = 0.5) +
        geom_smooth() +
        facet_wrap(~variable, scales = "free") +
        theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
        labs(title = plotTitles
             , subtitle = paste0("Numeric variables plotted against response variable (",Y,")")
             )
      
    }
    
    res$pairs <- ggpairs(datExp %>%
                           dplyr::select(1,where(~n_distinct(.) < 15))
                         ) +
      theme(axis.text.x=element_text(angle=90, vjust=0.5))
    
    return(res)
    
  }
  
  rr <- function(Taxa,Common,rrDf) {
    
    print(Taxa)
    
    outFile <- fs::path(outDir,paste0("reporting-rate_",Taxa,".rds"))
    
    res <- list()
  
    geos <- length(unique(rrDf$geo2))
    
    plotTitles <- bquote(~italic(.(Taxa))*":" ~ .(Common))
    
    # # GLM
    # if(geos > 1) {
    #   
    #   res$modGLM <- stan_glmer(cbind(success,trials-success) ~ year*geo2 + (1|cell)
    #                   , data = rrDf
    #                   , family = binomial()
    #                   , chains = if(testing) testChains else useChains
    #                   , iter = if(testing) testIter else useIter
    #                   )
    #   
    # } else {
    #   
    #   res$modGLM <- stan_glmer(cbind(success,trials-success) ~ year + (1|cell)
    #                     , data = rrDf
    #                     , family = binomial()
    #                     , chains = if(testing) testChains else useChains
    #                     , iter = if(testing) testIter else useIter
    #                     )
    #   
    # }
    
    # GAM
    if(geos > 1) {
      
      res$mod <- stan_gamm4(cbind(success,trials - success) ~ s(year, k = 4) + geo2 + s(year, k = 4, by = geo2)
                            , data = rrDf
                            , family = binomial()
                            , random = ~(1|cell)
                            , chains = if(testing) testChains else useChains
                            , iter = if(testing) testIter else useIter
                            )
      
    } else {
      
      res$mod <- stan_gamm4(cbind(success,trials-success) ~ s(year, k = 4)
                            , data = rrDf
                            , random = ~ (1|cell)
                            , family = binomial()
                            , chains = if(testing) testChains else useChains
                            , iter = if(testing) testIter else useIter
                            )
      
    }
    
    
    res$pred <- mod$data %>%
      distinct(year,geo2) %>%
      dplyr::mutate(col = row.names(.)
                    , success = 0
                    , trials = 100 #round(median(res$modGAM$data$trials),0)
                    ) %>%
      dplyr::left_join(as_tibble(posterior_predict(mod
                                                   , newdata = .
                                                   , re.form = NA#insight::find_formula(res$modGAM)$random[[1]]
                                                   )
                                 ) %>%
                         tibble::rownames_to_column(var = "row") %>%
                         tidyr::gather(col,value,2:ncol(.))
                       ) %>%
      dplyr::mutate(rawValue = as.numeric(value)
                    , value = rawValue/trials
                    )
    
    res$res <- res$pred %>%
      dplyr::group_by(year,geo2) %>%
      dplyr::summarise(modMean = mean(value)
                       , modMedian = quantile(value, 0.5)
                       , modci90lo = quantile(value, 0.05)
                       , modci90up = quantile(value, 0.95)
                       , text = paste0(round(modMedian,2)," (",round(modci90lo,2)," to ",round(modci90up,2),")")
                       ) %>%
      dplyr::ungroup()
    
    res$plot <- res$mod$data %>%
      distinct(year,geo2) %>%
      tidybayes::add_fitted_draws(modGAM
                                  , n = 100
                                  , re_formula = NA
                                  ) %>%
      ggplot(aes(x = year)) +
        geom_line(aes(y = .value, group = .draw), alpha = 0.05) +
        geom_line(data = res
                  , aes(y = modMean)
                  , size = 2
                  ) +
        geom_point(data = rrDf %>%
                     dplyr::group_by(year,geo2) %>%
                     dplyr::summarise(success = sum(success)
                                      , trials = n()
                                      ) %>%
                     dplyr::ungroup() %>%
                     dplyr::mutate(prop = success/trials)
                   , aes(y = prop
                         , colour = trials
                         )
                   ) +
        facet_wrap(~geo2) +
        scale_colour_viridis_c() +
        labs(title = plotTitles
             , subtitle = "Thick line is median credible value for that year"
             )
    
    write_rds(res,outFile)
    
  }
  
  
  ll <- function(Taxa,Common,rrDf) {
    
    print(Taxa)
    
    outFile <- fs::path(outDir,paste0("list-length_",Taxa,".rds"))
    
    res <- list()
    
    geos <- length(unique(rrDf$geo2))
    
    plotTitles <- bquote(~italic(.(Taxa))*":" ~ .(Common))
    
    # GLM
    # if(geos > 1) {
    #   
    #   res$mod <- stan_glmer(cbind(success,trials-success) ~ year*geo2*log(listLength) + (1|cell)
    #                       , data = rrDf
    #                       , family = binomial()
    #                       , chains = if(testing) testChains else useChains
    #                       , iter = if(testing) testIter else useIter
    #                       )
    #   
    # } else {
    #   
    #   res$mod <- stan_glmer(cbind(success,trials-success) ~ year*log(listLength) + (1|cell)
    #                       , data = rrDf
    #                       , family = binomial()
    #                       , chains = if(testing) testChains else useChains
    #                       , iter = if(testing) testIter else useIter
    #                       )
    #   
    # }
    
    # GAM
    if(geos > 1) {
      
      res$mod <- stan_gamm4(cbind(success,trials - success) ~
                              s(year, k = 4) +
                              geo2 +
                              logListLength +
                              s(year, k = 4, by = geo2) +
                              s(year, k = 4, by = logListLength) +
                              geo2*logListLength
                            , data = rrDf %>%
                              dplyr::mutate(logListLength = log(listLength))
                            , family = binomial()
                            , random = ~(1|cell)
                            , chains = if(testing) testChains else useChains
                            , iter = if(testing) testIter else useIter
                            )
      
    } else {
      
      res$mod <- stan_gamm4(cbind(success,trials-success) ~ s(year, k = 4)
                            , data = rrDf
                            , random = ~ (1|cell)
                            , family = binomial()
                            , chains = if(testing) testChains else useChains
                            , iter = if(testing) testIter else useIter
      )
      
    }
    
    res$pred <- rrDf %>%
      dplyr::distinct(geo1,geo2,year) %>%
      dplyr::mutate(listLength = median(rrDf$listLength)
                    , logListLength = log(listLength)
                    , col = row.names(.)
                    , success = 0
                    , trials = 100
                    ) %>%
      dplyr::left_join(as_tibble(posterior_predict(res$mod
                                                   , newdata = .
                                                   , re.form = NA#insight::find_formula(res$mod)$random
                                                   )
                                 ) %>%
                         tibble::rownames_to_column(var = "row") %>%
                         tidyr::gather(col,value,2:ncol(.))
                       ) %>%
      dplyr::mutate(rawValue = as.numeric(value)
                    , value = rawValue/trials
                    )
    
    res$res <- res$pred %>%
      dplyr::group_by(year,geo1,geo2,listLength) %>%
      dplyr::summarise(n = n()
                       , nCheck = nrow(as_tibble(res$mod))
                       , modMean = mean(value)
                       , modMedian = quantile(value, 0.5)
                       , modci90lo = quantile(value, 0.05)
                       , modci90up = quantile(value, 0.95)
                       , text = paste0(round(modMedian,2)," (",round(modci90lo,2)," to ",round(modci90up,2),")")
                       ) %>%
      dplyr::ungroup()
    
    plotData <- rrDf %>%
      dplyr::distinct(geo1,geo2,year) %>%
      dplyr::mutate(success = 0
                    , trials = 100
                    ) %>%
      dplyr::full_join(tibble(probs = quantProbs[2]) %>%
                         dplyr::mutate(listLength = map_dbl(probs
                                                            ,~quantile(unique(rrDf$listLength)
                                                                       ,probs = .
                                                                       )
                                                            )
                                       , logListLength = log(listLength)
                                       )
                       , by = character()
                       ) %>%
      dplyr::mutate(length = paste0("At list length quantile ",probs," = ",listLength)) %>%
      tidybayes::add_fitted_draws(res$mod
                                  , n = 500
                                  , re_formula = NA
                                  )
    
    res$plotLine <- plotData %>%
      ggplot(aes(x = year, y = prop)) +
        geom_line(aes(y = .value, group = .draw), alpha = 0.05) +
        geom_jitter(data = rrDf %>%
                      dplyr::group_by(year,geo2,listLength) %>%
                      dplyr::summarise(success = sum(success)
                                       , trials = n()
                                       ) %>%
                      dplyr::ungroup() %>%
                      dplyr::mutate(prop = success/trials)
                    , aes(colour = listLength)
                    , height = 0.025
                    , width = 0.2
                    ) +
        facet_wrap(~geo2) +
        scale_colour_viridis_c() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        labs(title = plotTitles
             , subtitle = unique(plotData$length)
             )
    
    res$plotRibbon <- ggplot() +
      geom_ribbon(data = res$res
                  , aes(year,modMean,ymin = modci90lo, ymax = modci90up)
                  , alpha = 0.4
                  ) +
      geom_line(data = res$res
                , aes(year,modMean)
                , linetype = 1
                , size = 1.5
                ) +
      geom_point(data = rrDf %>%
                   dplyr::group_by(year,geo2,listLength) %>%
                   dplyr::summarise(success = sum(success)
                                    , trials = n()
                                    , prop = success/trials
                                    )
                 ,aes(year
                      ,prop
                      ,colour = if(sum(grepl("listlength",tolower(names(res$mod$data))))>0) listLength else "none"
                      )
                 ) +
      facet_wrap(~geo2) +
      scale_colour_viridis_c() +
      labs(title = plotTitles
           , subtitle = unique(plotData$length)
           , colour = "List length"
           )
    
    write_rds(res,outFile)
    
  }
  
  
  
  
  add_residuals <- function(mod) {
    
    res <- list()
    
    #plotTitles <- bquote(~italic(.(Taxa))*":" ~ .(Common))
    
    res$resid <- tibble(residual = residuals(mod)
           , fitted = fitted(mod)
           ) %>%
      dplyr::bind_cols(as_tibble(mod$data))
    
    res$residPlot <- ggplot(res$resid, aes(fitted,residual)) +
      geom_point(size = 2) +
      geom_hline(aes(yintercept = 0), linetype = 2, colour = "red") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) %>%
      scale_colour_viridis_d(end=0.9)
    
    hasNumeric <- res$resid %>%
      dplyr::select(-1) %>%
      dplyr::select(where(is.numeric)) %>%
      ncol() %>%
      `>` (0)
    
    hasCharacter <- res$resid %>%
      dplyr::select(-1) %>%
      dplyr::mutate(across(where(is.factor),as.character)) %>%
      dplyr::select(where(is.character)) %>%
      ncol() %>%
      `>` (0)
    
    
    res$residPlotNum <- if(hasNumeric) {
      
      ggplot(res$resid %>%
               dplyr::select_if(is.numeric) %>%
               tidyr::pivot_longer(2:ncol(.))
             , aes(value,residual)
             ) +
        geom_point(size = 2) +
        geom_smooth(method = "lm")  +
        geom_hline(aes(yintercept = 0), linetype = 2, colour = "red") +
        facet_wrap(~name
                   , scales = "free_x"
                   ) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_colour_viridis_d()
      
    } else NULL
    
    res$residPlotChar <- if(hasCharacter) {
      
      ggplot(res$resid %>%
               dplyr::mutate(across(where(is.factor),as.character)) %>%
               dplyr::select(1,where(is.character)) %>%
               tidyr::pivot_longer(2:ncol(.))
             , aes(value,residual)
             ) +
        geom_boxplot() +
        geom_hline(aes(yintercept = 0), linetype = 2, colour = "red") +
        facet_wrap(~name, scales = "free") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      
    } else NULL
    
    return(res)
    
  }
  
  make_year_effect_df <- function(mod) {
    
   geos <- length(unique(mod$data$geo2))
   
   geosNames <- unique(mod$data$geo2)
   
   if("listLength" %in% names(mod$data)) atLevel <- median(mod$data$listLength)
   
   reference <- setdiff(unique(mod$data$geo2)
                        , str_extract(names(as_tibble(mod)),paste0(geosNames,collapse="|"))
                        )
   
   resRef <- as_tibble(mod) %>%
     dplyr::select(contains("year")) %>%
     dplyr::select(grep(paste0(geosNames,collapse="|"),names(.),value=TRUE,invert = TRUE)) %>%
     dplyr::rowwise() %>%
     dplyr::mutate(geo2 = reference
                   , across(grep("log",names(.)),~log(atLevel)*.)
                   , yearEff = sum(across(where(is.double)))
                   ) %>%
     dplyr::ungroup() %>%
     tibble::rowid_to_column()
   
   if(geos > 1) {
     
     resNotRef <- as_tibble(mod) %>%
       dplyr::select(contains("year")) %>%
       dplyr::rowwise() %>%
       dplyr::mutate(across(grep("log",names(.))
                            ,~log(atLevel)*.
                            )
                     , generic = sum(c_across(grep(paste0(geosNames,collapse = "|"),names(.),value=TRUE,invert = TRUE)))
                     ) %>%
       dplyr::ungroup() %>%
       dplyr::select(any_of(c("generic",grep(paste0(geosNames,collapse="|"),names(.),value=TRUE)))) %>%
       tidyr::pivot_longer(2:ncol(.)) %>%
       dplyr::mutate(geo2 = str_extract(name,paste0(geosNames,collapse="|"))
                     , term = gsub(paste0(geosNames,collapse="|"),"",name)
                     ) %>%
       dplyr::select(-name) %>%
       tidyr::pivot_wider(names_from = "term", values_from = "value") %>%
       tibble::rowid_to_column() %>%
       dplyr::rowwise() %>%
       dplyr::mutate(yearEff = sum(c_across(where(is.double)))) %>%
       dplyr::ungroup()
     
   }
   
   res <- resRef %>%
     purrr::when(geos > 1 ~ (.) %>% dplyr::bind_rows(resNotRef)
                 , geos == 1 ~ (.)
                 ) %>%
     dplyr::select(rowid,geo2,yearEff)
   
   return(res)
  
  }
    
  year_effect <- function(yearEffectDf,groups = "geo2") {
    
    nCheck <- yearEffectDf %>%
      dplyr::group_by(!!ensym(groups)) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::pull(n) %>%
      unique() %>%
      mean() # hack to deal with odd situation with, say, 3999 estimates
    
    yearEffectDf %>%
      dplyr::select(any_of(c("yearEff",groups))) %>%
      dplyr::group_by(!!ensym(groups)) %>%
      dplyr::summarise(n = n()
                       , nCheck = nCheck
                       , increasing = sum(yearEff > 0)/nCheck
                       , decreasing = sum(yearEff < 0)/nCheck
                       , meanEff = mean(yearEff)
                       , medianEff = median(yearEff)
                       , cilo = quantile(yearEff, probs = 0.05)
                       , ciup = quantile(yearEff, probs = 0.95)
                       ) %>%
      dplyr::mutate(likelihood = map(decreasing
                                     , ~cut(.
                                            , breaks = c(0,luLikelihood$maxVal)
                                            , labels = luLikelihood$likelihood
                                            , include.lowest = TRUE
                                            )
                                     )
                    ) %>%
      tidyr::unnest(cols = c(likelihood)) %>%
      dplyr::mutate(text = paste0(tolower(likelihood)
                                  , " to be decreasing in the "
                                  , geo2
                                  , " IBRA Subregion ("
                                  , 100*round(decreasing,2)
                                  , "% chance)"
                                  )
                    )
    
  }
  
  year_effect_plot <- function(yearEffectDf) {
    
    yearEffectDf %>%
      dplyr::group_by(geo2) %>%
      dplyr::mutate(decreasing = sum(yearEff<0)/n()) %>%
      ggplot(aes(yearEff,geo2, fill = decreasing)) +
        geom_density_ridges() +
        geom_vline(aes(xintercept = 0)
                   , linetype = 2
                   , colour = "red"
                   ) +
        scale_fill_viridis_c(limits = c(0,1)) +
        labs(x = "Effect of year"
             , y = "IBRA Subregion"
             , fill = "Probability of decline"
             )
    
  }
  
  make_year_comparison_df <- function(mod,years = c(2000,2020)) {
    
    res$mod$data %>%
      dplyr::distinct(geo1,geo2) %>%
      dplyr::full_join(tibble(year = years)
                       , by = character()
                       ) %>%
      dplyr::mutate(listLength = median(res$mod$data$listLength)
                    , col = row.names(.)
                    ) %>%
      dplyr::left_join(as_tibble(posterior_predict(res$mod
                                                   , newdata = .
                                                   , re.form = NA#insight::find_formula(res$mod)$random
                                                   , type = "response"
                                                   )
                                 ) %>%
                         tibble::rownames_to_column(var = "row") %>%
                         tidyr::gather(col,value,2:ncol(.))
                       )
    
    
  }
  
  overall_year_effect <- function(dfRow) {
    
    dfRow %>%
      tidyr::pivot_longer(1:ncol(.)) %>%
      tidyr::unnest(cols = c(value)) %>%
      dplyr::summarise(n = n()
                       , increase = sum(yearEff > 0)/n
                       , decline = sum(yearEff < 0)/n
                       , meanEff = mean(yearEff)
                       , medianEff = median(yearEff)
                       , cilo = quantile(yearEff, probs = 0.05)
                       , ciup = quantile(yearEff, probs = 0.95)
                       ) %>%
      dplyr::mutate(likelihood = map(decline
                                     , ~cut(.
                                            , breaks = c(0,luLikelihood$maxVal)
                                            , labels = luLikelihood$likelihood
                                            , include.lowest = TRUE
                                            )
                                     )
                    ) %>%
      tidyr::unnest(cols = c(likelihood)) %>%
      dplyr::mutate(text = paste0(tolower(likelihood)
                                  , " to be declining ("
                                  , 100*round(decline,2)
                                  , "% chance)"
                                  )
                    )
    
  }
 
 #------Run models-------
 
  
  # Check if rr models have been run - run if not
  todo <- datRR %>%
    dplyr::mutate(outFile = fs::path(outDir,paste0("reporting-rate_",Taxa,".rds"))
                  , done = map_lgl(outFile,file.exists)
                  ) %>%
    dplyr::filter(!done)
  
  if(nrow(todo) > 0) {
    
    future_pwalk(list(todo$Taxa
                      ,todo$Common
                      ,todo$data
                      )
               ,rr
               )
    
  }
  
  
  # Check if ll models have been run - run if not
  todo <- datRR %>%
    dplyr::mutate(outFile = fs::path(outDir,paste0("list-length_",Taxa,".rds"))
                  , done = map_lgl(outFile,file.exists)
                  ) %>%
    dplyr::filter(!done)
  
  if(nrow(todo) > 0) {
    
    future_pwalk(list(todo$Taxa
                        ,todo$Common
                        ,todo$data
                        )
                   ,ll
                 )
    
  }
  
  
  # Import results of models
  taxaModsFull <- datRR %>%
    dplyr::mutate(rr = fs::path(outDir,paste0("reporting-rate_",Taxa,".rds"))
                  , ll = fs::path(outDir,paste0("list-length_",Taxa,".rds"))
                  , exists = map2_lgl(rr,ll,~(file.exists(.x)&file.exists(.y)))
                  ) %>%
    dplyr::filter(exists) %>%
    dplyr::mutate(rr = map(rr,read_rds)
                  , ll = map(ll,read_rds)
                  , rrExp = pmap(list(Taxa
                                      , Common
                                      , data
                                      )
                                 ,explore
                                 )
                  # , llExp = pmap(list(Taxa
                  #                     , Common
                  #                     , dataLL
                  #                     ,"success"
                  #                     ,30
                  #                     )
                  #                ,explore
                  #                )
                  , rrResid = map(rr,~add_residuals(.$mod))
                  , llResid = map(ll,~add_residuals(.$mod))
                  , rrYearEffDf = map(rr,~.$mod %>% make_year_effect_df)
                  , llYearEffDf = map(ll,~.$mod %>% make_year_effect_df)
                  , rrYearEff = map(rrYearEffDf,year_effect,groups = "geo2")
                  , llYearEff = map(llYearEffDf,year_effect,groups = "geo2")
                  , rrYearEffPlot = map(rrYearEffDf,year_effect_plot)
                  , llYearEffPlot = map(llYearEffDf,year_effect_plot)
                  )
   
  taxaModsFull$overall <- taxaModsFull %>%
   dplyr::select(contains("YearEffDf")) %>%
   split(seq(nrow(.))) %>%
   purrr::map(overall_year_effect) %>%
   tibble() %>%
   tidyr::unnest(cols = 1) %>%
   dplyr::pull(text)

  
  timer$stop("rrll", comment = paste0("models run for ",nrow(taxaModsFull)," taxa"))
  