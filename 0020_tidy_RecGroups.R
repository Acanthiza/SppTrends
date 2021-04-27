
  timer$start("recGroups")
  
  #-------Filter for analysis---------
  
  
  
  datGroups <- datTidy %>%
    tidyr::nest(data = -c(!!ensym(taxGroup),geo1,geo2)) %>%
    dplyr::mutate(check = map_dbl(data,nrow)) %>%
    dplyr::filter(check > 2) %>%
    dplyr::mutate(check = map_dbl(data,~n_distinct(.$Taxa))) %>%
    dplyr::filter(check > 2*min(possibleGroups)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(groups = list(make_rec_groups(!!ensym(taxGroup),geo2,data)))
    
  
    
  
 
   
  timer$stop("recGroups")
  