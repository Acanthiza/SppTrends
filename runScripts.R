
  outName <- "2021-03-15-1829_the agricultural zone of South Australia"
    
  library(magrittr)
  
  runFrom <- 0
  runTo <- 50
  excludes <- 0
  
  testing <- T
  testRmd <- F
  commitToGit <- T
  
  getNewData <- FALSE
  
  runReason <- paste0("Slightly broader testing of RR, LL and OC. Tweaked filter_taxa_data. Testing = ",testing,".")
  
  dir() %>%
    grep("^\\d{4}_.*\\.R$",.,value=TRUE) %>%
    #grep(excludes,.,value = TRUE,invert = TRUE) %>%
    setNames(stringr::str_extract(.,"\\d{4}")) %>%
    `[` (names(.)[as.numeric(names(.)) <= runTo & as.numeric(names(.)) >= runFrom & !as.numeric(names(.)) %in% excludes]) %>%
    purrr::walk(source, verbose = TRUE)
  