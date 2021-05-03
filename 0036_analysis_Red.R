
  
  timer$start("red")
  
  #-------Filter for analysis---------
  
  red <- datTidy %>%
    tidyr::nest(data = -c(Taxa,!!ensym(taxGroup))) %>%
    #dplyr::sample_n(3) %>% # TESTING
    dplyr::mutate(datMat = map(data,~as.matrix(cbind(.$LONGITUDE,.$LATITUDE)))
                  , aoo = map_dbl(datMat,red::aoo)
                  , eoo = map_dbl(datMat,red::eoo)
                  )
    
  timer$stop("red",comment = paste0("Red list method used to generate area of occupancy and extent of occurence for "
                                    , nrow(red)
                                    , " taxa"
                                    )
             )
