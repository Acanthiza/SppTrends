

  timer$start("occupancy")
  
  #-------datForOcc---------
  
  # Sampling unit for occupancy is based on a cell within a yearmon. Thus, success is presence within a cell in a yearmon
  
  datForOcc <- dat %>%
    dplyr::mutate(quart = cut(month,breaks = c(0,3,6,9,12),labels = FALSE)
                  , yearquart = paste0(year,quart)
                  ) %>%
    dplyr::distinct(Taxa,year,quart,yearquart,!!ensym(taxGroup),geo1,geo2,cell) %>%
    dplyr::mutate(list = paste0(yearquart,"-",!!ensym(taxGroup),"-",cell)) %>%
    dplyr::mutate(success = 1) %>%
    dplyr::add_count(list, name = "listLength") %>%
    filter_taxa_data(minListLengthThresh = 1)
  
  taxaGeo <- datForOcc %>%
    dplyr::distinct(!!ensym(taxGroup)
                    ,Taxa
                    ,cell
                    ) 
  
  #-------datForRR---------
  
  datOcc <- datForOcc %>%
    dplyr::inner_join(taxaGeo) %>%
    dplyr::mutate(success = 1) %>%
    tidyr::nest(data = c(!!ensym(taxGroup),Taxa,cell,yearquart,success))
    tidyr::pivot_wider(names_from = "Taxa", values_from = "success", values_fill = 0) %>%
    tidyr::pivot_longer(cols = any_of(unique(taxaGeo$Taxa)),names_to = "Taxa", values_to = "success") %>%
    dplyr::inner_join(taxaGeo) %>%
    purrr::when(testing ~ (.) %>% dplyr::filter(Taxa %in% tests)
                , !testing ~ (.)
                ) %>%
    tidyr::nest(data = c(geo1,geo2,cell,year,quart,yearquart,list,listLength,success))
  
    %>%
      dplyr::left_join(luTax %>%
                         dplyr::select(!!ensym(taxGroup),Taxa,Common)
      )
  
  occ <- function(Taxa,Common,geo = "geo2",occDf) {
    
    print(Taxa)
    
    outFile <- fs::path(outDir,paste0("occupancy_",Taxa,".rds"))
    
    res <- list()
    
    occDf %>%
      tidyr::nest(data = -c(!!ensym(geo),year))
      
    
    
    
  }
  
  
  
  