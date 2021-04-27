  
  
# What other taxa are recorded in conjunction with taxa X? These help create absences....

  make_clusters <- function(data
                            ,methodsDf = clustMethod
                            ,sppCol = "Taxa"
                            ,siteCol = "list"
                            ,groups = possibleGroups[possibleGroups < 0.5*length(unique(data$Taxa))]
                            ,minTaxaCount = 1
                            ) {
    
    datWide <- data %>%
      dplyr::add_count(!!ensym(sppCol)) %>%
      dplyr::filter(n > minTaxaCount) %>%
      dplyr::select(!!ensym(sppCol),!!ensym(siteCol)) %>%
      dplyr::mutate(p = 1) %>%
      tidyr::pivot_wider(names_from = all_of(sppCol), values_from = "p", values_fill = 0)
    
    siteNames <- datWide %>% dplyr::pull(!!ensym(siteCol))
    
    dist <- parDist(datWide %>% tibble::column_to_rownames(siteCol) %>% as.matrix()
                    , method = "bray"
                    , threads = useCores
                    )
    
    dend <- methodsDf %>%
      dplyr::mutate(dend = map(method
                               ,~fastcluster::hclust(dist, .)
                               )
                    )
    
    clust <- dend %>%
      dplyr::mutate(clusters = map(dend
                                   , cutree
                                   , groups
                                   )
                    , clusters = map(clusters
                                     , as_tibble
                                     )
                    ) %>%
      dplyr::select(-dend) %>%
      tidyr::unnest(clusters) %>%
      tidyr::pivot_longer(2:ncol(.),names_to = "groups",values_to ="clust") %>%
      dplyr::mutate(groups = as.integer(groups)) %>%
      tidyr::nest(clusters = c(clust)) %>%
      dplyr::mutate(clusters = future_map(clusters
                                          , . %>%
                                            dplyr::mutate(!!ensym(siteCol) := siteNames
                                                          , cluster = numbers2words(clust)
                                                          , cluster = fct_reorder(cluster,clust)
                                                          )
                                          )
                    )
    
  }
  
  clusters_summarise <- function(clustersDf) {
    
    clustersDf %>%
      dplyr::mutate(summaries = future_map(clusters,cluster_summarise)) %>%
      tidyr::unnest(cols = c(summaries)) %>%
      dplyr::mutate(sil = future_map(clusters,make_sil,distObj = dist)
                    , silDf = future_map2(clusters,sil,make_sil_df)
                    , macroSil = map_dbl(silDf,~mean(.$sil_width))
                    , indval = future_map(clusters,cluster_indval,datWide)
                    , indValDf = future_map2(indval,clusters,cluster_indval_df)
                    , wss = future_map2(clusters,clusters,calc_SS,dist = dist)
                    , macroWSS = map_dbl(wss,~sum(.$wss))
                    )
    
    clustDiagnostics <- diagnostic_df(clust, topThresh = 0.9, bestThresh = 1, diagnosticMinGroups = min(possibleGroups))
    
    # diagnostic_plot(clustDiagnostics)
    
    clustBest <- clustDiagnostics %>%
      dplyr::filter(best) %>%
      dplyr::distinct(method,groups) %>%
      dplyr::inner_join(clust %>%
                          dplyr::select(method,groups,indValDf,contains("macro"))
                        ) %>%
      tidyr::unnest(cols = c(indValDf)) %>%
      dplyr::filter(macroSil == max(macroSil)) %>%
      dplyr::filter(macroWSS == min(macroWSS))
    
  }



  add_list_length <- function(df,groups) {
    
    df %>%
      dplyr::mutate(list = paste(!!!rlang::syms(groups),sep="_")) %>%
      dplyr::add_count(list, name = "listLength") %>%
      dplyr::mutate(success = 1)
    
  }
  
    
  mod_explore <- function(Taxa
                          , Common
                          , df
                          , modPath
                          , modType
                          , respVar = "prop"
                          , expVar = alwaysGroup %>% grep("yday|site",.,value = TRUE,invert = TRUE)
                          , maxLevels = 30
                          , draws = 200
                          , postGroups = c("Taxa"
                                           ,"Common"
                                           ,"listLength"
                                           ,"year"
                                           ,"geo2"
                                           )
                          ) {
    
    print(Taxa)
    
    #-------import/export--------
    mod <- read_rds(modPath)
    
    outPred <- path(outDir,paste0(gsub(" ","-",tolower(modType)),"_Pred_",Taxa,".feather"))
    
    outRes <- path(outDir,paste0(gsub(" ","-",tolower(modType)),"_Res_",Taxa,".rds"))
    
    outYearDiff <- path(outDir,paste0(gsub(" ","-",tolower(modType)),"_YearDif_",Taxa,".feather"))
    
    #-------setup explore-------
    
    res <- list()
    
    plotTitles <- bquote(~italic(.(Taxa))*":" ~ .(Common))
    
    hasLL <- sum(grepl("listLength",names(mod$coefficients))) > 0
    
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
                       , stat = "count") +
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
                             tidyr::gather(variable,value,2:ncol(.)) %>%
                             dplyr::group_by(variable) %>%
                             dplyr::mutate(levels = n_distinct(value)) %>%
                             dplyr::ungroup() %>%
                             dplyr::filter(levels < maxLevels)
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
                           dplyr::mutate(across(where(is.character),factor)) %>%
                           dplyr::select(where(~is.numeric(.x)|is.factor(.x) & n_distinct(.x) < 15)) %>%
                           dplyr::mutate(across(where(is.factor),factor))
                         ) +
      theme(axis.text.x=element_text(angle=90, vjust=0.5))
    
    
    #-------residuals-------
    
    if(length(residuals(mod)) == nrow(df)) {
      
      res$resid <- tibble(residual = residuals(mod)
                        , fitted = fitted(mod)
                        ) %>%
      dplyr::bind_cols(df)
      

      res$residPlot <- ggplot(res$resid, aes(fitted,residual)) +
        geom_point(size = 2) +
        geom_hline(aes(yintercept = 0), linetype = 2, colour = "red") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) %>%
        scale_colour_viridis_d(end=0.9)
      
  
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
                 tidyr::pivot_longer(2:ncol(.)) %>%
                 dplyr::group_by(name) %>%
                 dplyr::mutate(levels = n_distinct(value)) %>%
                 dplyr::ungroup() %>%
                 dplyr::filter(levels < maxLevels)
               , aes(value,residual)
               ) +
          geom_boxplot() +
          geom_hline(aes(yintercept = 0), linetype = 2, colour = "red") +
          facet_wrap(~name, scales = "free") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

      } else NULL
      
    }
    


    #---------post explore-------
    
    if(family(mod)$family == "beta") class(mod) <- unique(c(class(mod),"betareg"))
    
    isBinomialMod <- family(mod)$family == "binomial"

    pred <- df %>%
      dplyr::distinct(geo2) %>%
      dplyr::left_join(df %>%
                         dplyr::distinct(year)
                       , by = character()
                       ) %>%
      dplyr::mutate(listLength = if(hasLL) median(df$listLength) else NULL
                    , listLengthLog = if(hasLL) log(listLength) else NULL
                    , col = row.names(.)
                    , success = if(isBinomialMod) 0 else NULL
                    , trials = if(isBinomialMod) 100 else NULL
                    ) %>%
      dplyr::left_join(as_tibble(posterior_predict(mod
                                                   , newdata = .
                                                   , re.form = NA#insight::find_formula(mod)$random
                                                   , type = "response"
                                                   )
                                 ) %>%
                         tibble::rownames_to_column(var = "row") %>%
                         tidyr::gather(col,value,2:ncol(.))
                       ) %>%
      dplyr::mutate(rawValue = as.numeric(value)
                    , value = if(isBinomialMod) rawValue/trials else rawValue
                    )
    
    write_feather(pred,outPred)

    res$res <- pred %>%
      dplyr::group_by(across(any_of(postGroups))) %>%
      dplyr::summarise(n = n()
                       , nCheck = nrow(as_tibble(mod))
                       , modMean = mean(value)
                       , modMedian = quantile(value, 0.5)
                       , modci90lo = quantile(value, 0.05)
                       , modci90up = quantile(value, 0.95)
                       , text = paste0(round(modMedian,2)," (",round(modci90lo,2)," to ",round(modci90up,2),")")
                       ) %>%
      dplyr::ungroup()


    #------res plot data-------

    plotData <- df %>%
      dplyr::distinct(geo2,year) %>%
      dplyr::mutate(success = 0
                    , trials = 100
                    ) %>%
      dplyr::full_join(tibble(probs = quantProbs[2]) %>%
                         {if(hasLL) (.) %>% dplyr::mutate(listLength = map_dbl(probs
                                                                               ,~quantile(unique(df$listLength)
                                                                                          ,probs = .
                                                                                          )
                                                                               )
                                                          , listLengthLog = log(listLength)
                                                          , length = paste0("At list length quantile ",probs," = ",listLength)
                         ) else (.)
                           }
                       , by = character()
                       ) %>%
      tidybayes::add_fitted_draws(mod
                                  , n = draws
                                  , re_formula = NA
                                  )

    subTitle <-  if(hasLL) {

      paste0("List length corrected reporting rate.\nDashed red lines indicate years for comparison (see text).")

    } else {

      paste0(modType,".\nDashed red lines indicate years for comparison (see text).")

    }

    subTitleLine <- paste0(
      subTitle
      , if(hasLL) {

        paste0("\nLines are ",draws," draws from posterior distribution.\n",unique(plotData$length))

      } else {

        paste0("\nLines are ",draws," draws from posterior distribution.")

      }
    )

    subTitleRibbon <- paste0(subTitle,"\nMedian (thick line) and 90% credible intervals (shaded).")

    #-------res plotLine-----------

    p <- plotData %>%
      ggplot(aes(x = year, y = !!ensym(respVar))) +
      geom_line(aes(y = .value, group = .draw), alpha = 0.5) +
      geom_vline(xintercept = testYears$year, linetype = 2, colour = "red") +
      facet_wrap(~geo2) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      labs(title = plotTitles
           , subtitle = subTitleLine
           )

    if(hasLL) p <- p +
      geom_jitter(data = df
                  ,aes(year
                       , !!ensym(respVar)
                       , colour = listLength
                       )
                  , width = 0.1
                  , height = 0.05
                  ) +
      scale_colour_viridis_c() +
      labs(colour = "List length")

    if(!hasLL) p <- p +
      geom_jitter(data = df
                  ,aes(year
                       , !!ensym(respVar)
                       , colour = trials
                       )
                  , width = 0.1
                  , height = 0.01
                  ) +
      scale_colour_viridis_c()

    res$plotLine <- p


    #------res plotRibbon-------

    p <- ggplot() +
      geom_ribbon(data = res$res
                  , aes(year,modMean,ymin = modci90lo, ymax = modci90up)
                  , alpha = 0.4
                  ) +
      geom_line(data = res$res
                , aes(year,modMean)
                , linetype = 1
                , size = 1.5
                ) +
      geom_vline(xintercept = testYears$year, linetype = 2, colour = "red") +
      facet_wrap(~geo2) +
      labs(title = plotTitles
           , subtitle = subTitleRibbon
           )

    if(hasLL) p <- p +
      geom_jitter(data = df
                  ,aes(year
                       ,!!ensym(respVar)
                       , colour = listLength
                       )
                  , width = 0.1
                  , height = 0.05
                  ) +
      scale_colour_viridis_c() +
      labs(colour = "List length")

    if(!hasLL) p <- p +
      geom_jitter(data = df
                  ,aes(year
                       , !!ensym(respVar)
                       , colour = trials
                       )
                  , width = 0.1
                  , height = 0.05
                  ) +
      scale_colour_viridis_c() +
      labs(colour = "Trials")

    res$plotRibbon <- p


    #------year difference df-----------

    res$yearDifferenceDf <- df %>%
      dplyr::distinct(across(any_of(postGroups[postGroups != "year" & postGroups != "listLength"]))) %>%
      dplyr::full_join(testYears
                       , by = character()
                       ) %>%
      dplyr::mutate(listLength = if(hasLL) median(df$listLength) else NULL
                    , listLengthLog = if(hasLL) log(listLength) else NULL
                    , col = row.names(.)
                    , success = 0
                    , trials = 100
                    , nCheck = nrow(as_tibble(mod))
                    , modType = modType
                    , Taxa = Taxa
                    , Common = Common
                    ) %>%
      dplyr::left_join(as_tibble(posterior_predict(mod
                                                   , newdata = .
                                                   , re.form = NA#insight::find_formula(mod)$random
                                                   )
                                 ) %>%
                         tibble::rownames_to_column(var = "row") %>%
                         tidyr::gather(col,value,2:ncol(.))
                       ) %>%
      dplyr::select(-c(col)) %>%
      dplyr::mutate(value = as.numeric(if(isBinomialMod) value/trials else value)) %>%
      tidyr::pivot_wider(names_from = "type"
                         , values_from = c("year","value")
                         ) %>%
      setNames(gsub("\\d{4}","",names(.))) %>%
      dplyr::mutate(diff = as.numeric(value_recent-value_reference))
    
    write_feather(res$yearDifferenceDf,outYearDiff)


    #-------year difference res---------

    res$yearDifferenceRes <- res$yearDifferenceDf %>%
      dplyr::group_by(across(any_of(c("type","modType",postGroups)))) %>%
      dplyr::summarise(n = n()
                       , nCheck = unique(nCheck)
                       , lower = sum(diff < 0)/nCheck
                       , higher = sum(diff > 0)/nCheck
                       , meanDiff = mean(diff)
                       , medianDiff = median(diff)
                       , cilo = quantile(diff, probs = 0.05)
                       , ciup = quantile(diff, probs = 0.95)
                       , reference = unique(year_reference)
                       , recent = unique(year_recent)
                       ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(likelihood = map(lower
                                     , ~cut(.
                                            , breaks = c(0,luLikelihood$maxVal)
                                            , labels = luLikelihood$likelihood
                                            , include.lowest = TRUE
                                            )
                                     )
                    ) %>%
      tidyr::unnest(cols = c(likelihood)) %>%
      dplyr::mutate(text = paste0(tolower(likelihood)
                                  , " to be lower in "
                                  , geo2
                                  , " IBRA Subregion ("
                                  , 100*round(lower,2)
                                  , "% chance)"
                                  )
                    , text = gsub("in Kangaroo Island","on Kangaroo Island",text)
                    )

    #------year difference plot--------

    res$yearDifferencePlot <- res$yearDifferenceDf %>%
      dplyr::group_by(geo2) %>%
      dplyr::mutate(lower = sum(diff<0)/nCheck) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(likelihood = map(lower
                                     , ~cut(.
                                            , breaks = c(0,luLikelihood$maxVal)
                                            , labels = luLikelihood$likelihood
                                            , include.lowest = TRUE
                                            )
                                     )
                    ) %>%
      tidyr::unnest(cols = c(likelihood)) %>%
      dplyr::mutate(likelihood = fct_expand(likelihood,levels(luLikelihood$likelihood))) %>%
      ggplot(aes(diff,geo2,fill = likelihood)) +
      geom_density_ridges() +
      geom_vline(aes(xintercept = 0)
                 , linetype = 2
                 , colour = "red"
                 ) +
      scale_fill_viridis_d(drop = FALSE) +
      labs(title = plotTitles
           , subtitle = paste0("Difference in "
                               ,recent
                               ," "
                               ,tolower(modType)
                               ," compared to "
                               ,reference
                               )
           , x = "Difference"
           , y = "IBRA Subregion"
           , fill = "Likelihood of decrease"
           , caption = paste0("Red dotted line indicates no change from ",reference)
           )
    
    
    write_rds(res,outRes)
    
  }
  
  year_difference_overall <- function(Taxa,Common,yearDiffDfs) {
    
    plotTitles <- bquote(~italic(.(Taxa))*":" ~ .(Common))
    
    res <- list()
    
    res$yearDiffOverallRes <- yearDiffDfs %>%
      dplyr::summarise(n = n()
                       , increase = sum(diff > 0)/n
                       , decline = sum(diff < 0)/n
                       , meanEff = mean(diff)
                       , medianEff = median(diff)
                       , cilo = quantile(diff, probs = 0.05)
                       , ciup = quantile(diff, probs = 0.95)
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
    
    res$yearDiffText <- paste0(
      Taxa
      , if(!is.null(Common)) paste0(" (",Common,")")
      , " was "
      , res$yearDiffOverallRes$text
      , " across "
      , aoiFullName
      , " based on "
      , n_distinct(yearDiffDfs$type)
      , " models ("
      , vec_to_sentence(unique(yearDiffDfs$type))
      , ") using data from "
      , n_distinct(yearDiffDfs$geo2)
      , " IBRA Subregions ("
      , vec_to_sentence(unique(yearDiffDfs$geo2))
      , ")"
    )
    
    res$yearDiffOverallPlot <- yearDiffDfs %>%
      dplyr::mutate(likelihood = res$yearDiffOverallRes$likelihood) %>%
      ggplot(aes(diff,fill = likelihood)) +
      geom_density() +
      geom_vline(aes(xintercept = 0)
                 , linetype = 2
                 , colour = "red"
                 ) +
      scale_fill_viridis_d(drop = FALSE) +
      labs(title = plotTitles
           , subtitle = paste0("Distribution of credible values for change between "
                               ,recent
                               ," and "
                               ,reference
                               )
           , x = "Difference"
           , y = "IBRA Subregion"
           , fill = "Likelihood of decrease"
           , caption = paste0("Red dotted line indicates no change from ",reference)
      )
    
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
  
  

  
  filter_taxa_data <- function(df
                               , minListLengthThresh = 3
                               , maxListLengthOccurenceThresh = 3
                               , minlistOccurenceThresh = 5
                               , minYearsThresh = 3
                               , minListLengthsThresh = 2
                               , minYearSpanThresh = 20
                               ) {
    
    find_min_list_length <- function(df) {
      
      min(df$listLength)
      
    }
    
    find_max_list_length_occurence <- function(df) {
      
      df %>%
        dplyr::count(geo2,Taxa,listLength,name="lists") %>%
        dplyr::group_by(Taxa) %>%
        dplyr::filter(listLength == max(listLength)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(listLength == min(listLength)) %>%
        dplyr::pull(listLength) %>%
        unique()
      
    }
    
    find_min_list_occurence <- function(df) {
      
      df %>%
        dplyr::count(geo2,!!ensym(taxGroup),Taxa,name = "lists") %>%
        dplyr::filter(lists == min(lists)) %>%
        dplyr::pull(lists) %>%
        unique()
      
    }
    
    find_min_years <- function(df) {
      
      df %>%
        dplyr::count(geo2,!!ensym(taxGroup),Taxa,year,name = "blah") %>%
        dplyr::count(geo2,!!ensym(taxGroup),Taxa, name = "years") %>%
        dplyr::filter(years == min(years)) %>%
        dplyr::pull(years) %>%
        unique()
      
    }
    
    find_min_list_lengths <- function(df) {
      
      df %>%
        dplyr::count(geo2,Taxa,listLength,name = "blah") %>%
        dplyr::count(geo2,Taxa, name = "lengths") %>%
        dplyr::filter(lengths == min(lengths)) %>%
        dplyr::pull(lengths) %>%
        unique()
      
    }
    
    find_min_year_span <-function(df) {
      
      df %>%
        dplyr::group_by(!!ensym(taxGroup),geo2,Taxa) %>%
        dplyr::summarise(minYear = min(year)
                      , maxYear = max(year)
                      , diffYear = maxYear - minYear
                      ) %>%
        dplyr::ungroup() %>%
        dplyr::pull(diffYear) %>%
        min() %>%
        unique()
      
    }
    
    minListLength <- find_min_list_length(df)
    maxListLengthOccurence <- find_max_list_length_occurence(df)
    minlistOccurence <- find_min_list_occurence(df)
    minYears <- find_min_years(df)
    minLengths <- find_min_list_lengths(df)
    minYearSpan <- find_min_year_span(df)
    
    while(minListLength < minListLengthThresh |
          maxListLengthOccurence < maxListLengthOccurenceThresh |
          minlistOccurence < minlistOccurenceThresh |
          minYears < minYearsThresh |
          minLengths < minListLengthsThresh |
          minYearSpan < minYearSpanThresh
    ) {
      
      df <- df %>%
        dplyr::filter(listLength > minListLengthThresh)
      
      removeTaxaOnShortLists <- df %>%
        dplyr::count(geo2,!!ensym(taxGroup),Taxa,listLength,name="lists") %>%
        dplyr::group_by(Taxa) %>%
        dplyr::filter(listLength == max(listLength)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(listLength < maxListLengthOccurenceThresh) %>%
        dplyr::distinct(geo2,!!ensym(taxGroup),Taxa)
      
      removeTaxaWithFewOccurrences <- df %>%
        dplyr::anti_join(removeTaxaOnShortLists) %>%
        dplyr::count(geo2,!!ensym(taxGroup),Taxa,name = "lists") %>%
        dplyr::filter(lists < minlistOccurenceThresh) %>%
        dplyr::distinct(geo2,!!ensym(taxGroup),Taxa)
      
      removeTaxaWithFewYears <- df %>%
        dplyr::anti_join(removeTaxaOnShortLists) %>%
        dplyr::anti_join(removeTaxaWithFewOccurrences) %>%
        dplyr::count(geo2,!!ensym(taxGroup),Taxa,year,name = "blah") %>%
        dplyr::count(geo2,!!ensym(taxGroup),Taxa, name = "years") %>%
        dplyr::filter(years < minYearsThresh) %>%
        dplyr::distinct(geo2,!!ensym(taxGroup),Taxa)
      
      removeTaxaWithFewLengths <- df %>%
        dplyr::anti_join(removeTaxaOnShortLists) %>%
        dplyr::anti_join(removeTaxaWithFewOccurrences) %>%
        dplyr::anti_join(removeTaxaWithFewYears) %>%
        dplyr::count(geo2,Taxa,listLength,name = "blah") %>%
        dplyr::count(geo2,Taxa, name = "lengths") %>%
        dplyr::filter(lengths < minListLengthsThresh) %>%
        dplyr::distinct(geo2,Taxa)
      
      removeTooFewYears <- df %>%
        dplyr::anti_join(removeTaxaOnShortLists) %>%
        dplyr::anti_join(removeTaxaWithFewOccurrences) %>%
        dplyr::anti_join(removeTaxaWithFewYears) %>%
        dplyr::anti_join(removeTaxaWithFewLengths) %>%
        dplyr::group_by(!!ensym(taxGroup),geo2,Taxa) %>%
        dplyr::summarise(minYear = min(year)
                         , maxYear = max(year)
                         , diffYear = maxYear - minYear
                         ) %>%
        dplyr::ungroup() %>%
        dplyr::filter(diffYear < minYearSpanThresh) %>%
        dplyr::distinct(!!ensym(taxGroup),Taxa,geo2)
        
      
      df <- df %>%
        dplyr::anti_join(removeTaxaOnShortLists) %>%
        dplyr::anti_join(removeTaxaWithFewOccurrences) %>%
        dplyr::anti_join(removeTaxaWithFewYears) %>%
        dplyr::anti_join(removeTaxaWithFewLengths) %>%
        dplyr::anti_join(removeTooFewYears) %>%
        dplyr::add_count(list, name = "listLength")
      
      minListLength <- find_min_list_length(df)
      maxListLengthOccurence <- find_max_list_length_occurence(df)
      minlistOccurence <- find_min_list_occurence(df)
      minYears <- find_min_years(df)
      minLengths <- find_min_list_lengths(df)
      minYearSpan <- find_min_year_span(df)
      
      res <- paste0("Minimum list length = ",minListLength
                    ,"\nMaximum list taxa occurs on = ",maxListLengthOccurence
                    ,"\nMinimum list occurence = ",minlistOccurence
                    ,"\nMinimum years = ",minYears
                    ,"\nMinimum list lengths = ",minLengths
                    ,"\nMinimum year span = ",minYearSpan
                    ,"\nTotal records = ",nrow(df)
                    ,"\n"
                    )
      
      if(testing) cat(res)
      
    }
    
    return(df)
    
  }
  