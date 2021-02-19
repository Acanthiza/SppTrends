
  N <- 100
  
  df <- tibble(year = sample(1991:2020,N,replace = TRUE)) %>%
    dplyr::mutate(success = sample(c(0,1),N,replace = TRUE,prob = c(0.6,0.4))
                  , geo2 = sample(letters[1:5],N,replace = TRUE)
                  , cell = paste0(geo2,sample(1:10,N,replace = TRUE))
                  , listLength = rnbinom(100,5,mu = 15)
                  , listLength = if_else(listLength == 0, 1, listLength)
                  , logListLength = log(listLength)
                  , trials = 1
                  ) %>%
    dplyr::distinct()
  
  
  mod <- stan_gamm4(cbind(success,trials - success) ~
                      s(year, k = 4, bs = "ts") +
                      s(year, k = 4, by = geo2, bs = "ts") +
                      s(year, k = 4, by = logListLength, bs = "ts") +
                      geo2*logListLength
                    # + s(year, k = 4, by = c(geo2,logListLength), bs = "ts") # FAILS
                    
                    # https://stackoverflow.com/questions/47934100/how-to-specify-the-non-linear-interaction-of-two-factor-variables-in-generalised
                    # link used BAM, not GAM
                    # + s(year, k = 4, by = interaction(geo2,logListLength), bs = "ts") # FAILS 
                    , data = df
                    , family = binomial()
                    , random = ~(1|cell)
                    , chains = testChains
                    , iter = testIter
                    )
  
  
  pred <- mod$data %>%
    distinct(year,geo2) %>%
    dplyr::mutate(logListLength = log(median(rrDf$listLength))
                  , col = row.names(.)
                  , success = 0
                  , trials = 100 #round(median(res$mod$data$trials),0)
                  ) %>%
    dplyr::left_join(as_tibble(posterior_predict(mod
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
  
  res <- pred %>%
    dplyr::group_by(year,geo2) %>%
    dplyr::summarise(modMean = mean(value)
                     , modMedian = quantile(value, 0.5)
                     , modci90lo = quantile(value, 0.05)
                     , modci90up = quantile(value, 0.95)
                     , text = paste0(round(modMedian,2)," (",round(modci90lo,2)," to ",round(modci90up,2),")")
                     ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(medLL = exp(median(mod$data$logListLength)))
  
  mod$data %>%
    distinct(year,geo2) %>%
    dplyr::mutate(logListLength = log(median(rrDf$listLength))
                  , col = row.names(.)
                  , success = 0
                  , trials = 100 #round(median(res$mod$data$trials),0)
                  ) %>%
    tidybayes::add_fitted_draws(mod
                                , n = 100
                                , re_formula = NA
                                ) %>%
    ggplot(aes(x = year)) +
      geom_line(aes(y = .value, group = .draw), alpha = 0.05) +
      geom_line(data = res
                , aes(y = modMean)
                , size = 2
                ) +
      geom_point(data = df %>%
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
      scale_colour_viridis_c()
  
  ggplot() +
    geom_ribbon(data = res
                , aes(year,modMean,ymin = modci90lo, ymax = modci90up)
                , alpha = 0.4
                ) +
    geom_line(data = res
              , aes(year,modMean)
              , linetype = 1
              , size = 1.5
              ) +
    geom_point(data = df %>%
                 dplyr::group_by(year,geo2,listLength) %>%
                 dplyr::summarise(success = sum(success)
                                  , trials = n()
                                  , prop = success/trials
                                  )
               ,aes(year
                    ,prop
                    ,colour = if(sum(grepl("listlength",tolower(names(mod$data))))>0) listLength else "none"
                    )
               ) +
    facet_wrap(~geo2) +
    scale_colour_viridis_c() +
    labs(colour = "List length")
  