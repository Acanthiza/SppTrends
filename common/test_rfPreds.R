


  temp <- sampleRandom(rfPred,999,cells = TRUE) %>%
    as_tibble()
  
  tempVals <- cells_env(temp,unstack(rfPreds)) %>%
    tidyr::pivot_longer(2:ncol(.)) %>%
    dplyr::count(cell,value)
  