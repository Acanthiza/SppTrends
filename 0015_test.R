
if(testing) {
  
  #--------Test taxa---------
  
  # These are only used when testing = TRUE
  
  orders <- c(NULL
              # , "Anseriformes"
              # , "Charadriiformes"
              # , "Accipitriformes"
              # , "Diprotodontia"
              )
  
  genera <- c(NULL
              # , "Melithreptus"
              # , "Pandion"
              # , "Haliaeetus"
              # , "Cacatua"
              # , "Myiagra"
              # , "Pseudonaja"
              # , "Pogona"
              # , "Vespadelus"
              # , "Macropus"
              # , "Pachycephala"
              # , "Phylidonyris"
              # , "Climacteris"
              )
  
  
  kiSpp <- read_csv(path("data","kiSpp.csv")
                    , col_names = FALSE
                    ) %>%
    dplyr::mutate(X1 = gsub("\\s"," ",X1)) %>%
    tidyr::separate(X1,into = c("genus","species")) %>%
    dplyr::mutate(spp = paste0(genus," ",species)) %>%
    dplyr::distinct(spp) %>%
    dplyr::arrange(spp) %>%
    dplyr::pull(spp)
  
  
  spp <- c(NULL
          , "Melithreptus gularis"
          , "Melithreptus lunatus"
          )
  
  
  spp <- c(spp
           ,kiSpp
           )
  
    
  testOrder <- luTax %>%
    dplyr::filter(order %in% orders)
  
  testGenus <- luTax %>%
    dplyr::filter(genus %in% genera)
  
  testSpp <- luTax %>%
    dplyr::filter(Taxa %in% spp)
  
  tests <- sort(unique(c(testOrder$Taxa,testGenus$Taxa,testSpp$Taxa)))
  
  testTaxGroups <- luTax %>%
    dplyr::filter(Taxa %in% tests) %>%
    dplyr::pull(!!ensym(taxGroup)) %>%
    unique()
  
}

  
  
