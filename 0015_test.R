
if(testing) {
  
  #--------Test taxa---------
  
  # These are only used when testing = TRUE
  
  orders <- c(NULL
              #, "Anseriformes"
              #, "Charadriiformes"
              #, "Accipitriformes"
              )
  
  genera <- c(NULL
              #, "Melithreptus"
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
  
  testOrder <- luTax %>%
    dplyr::filter(Order %in% orders)
  
  testGenus <- luTax %>%
    dplyr::filter(Genus %in% genera)
  
  tests <- c(testOrder$Taxa,testGenus$Taxa)
  
}
  
