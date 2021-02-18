
if(testing) {
  
  #--------Test taxa---------
  
  # These are only used when testing = TRUE
  
  orders <- c(NULL
              #, "Anseriformes"
              #, "Charadriiformes"
              )
  
  genera <- c("Melithreptus"
              #, "Cacatua"
              , "Myiagra"
              #, "Pseudonaja"
              #, "Pogona"
              #, "Vespadelus"
              #, "Macropus"
              , "Pachycephala"
              #, "Phylidonyris"
              #, "Climacteris"
              )
  
  testOrder <- luTax %>%
    dplyr::filter(Order %in% orders)
  
  testGenus <- luTax %>%
    dplyr::filter(Genus %in% genera)
  
  tests <- c(testOrder$Taxa,testGenus$Taxa)
  
}
  
