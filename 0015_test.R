
if(testing) {
  
  #--------Test taxa---------
  
  # These are only used when testing = TRUE
  
  orders <- c(NULL
              #, "Anseriformes"
              #, "Charadriiformes"
              #, "Accipitriformes"
              #, "Diprotodontia"
              )
  
  genera <- c(NULL
              , "Melithreptus"
              , "Pandion"
              , "Haliaeetus"
              #, "Cacatua"
              #, "Myiagra"
              #, "Pseudonaja"
              #, "Pogona"
              #, "Vespadelus"
              #, "Macropus"
              #, "Pachycephala"
              #, "Phylidonyris"
              #, "Climacteris"
              )
  
  spp <- c(NULL
          , "Melithreptus gularis"
          , "Melithreptus lunatus"
          )
  
  testOrder <- luTax %>%
    dplyr::filter(Order %in% orders)
  
  testGenus <- luTax %>%
    dplyr::filter(Genus %in% genera)
  
  testSpp <- luTax %>%
    dplyr::filter(Taxa %in% spp)
  
  tests <- sort(unique(c(testOrder$Taxa,testGenus$Taxa,testSpp$Taxa)))
  
}
  
