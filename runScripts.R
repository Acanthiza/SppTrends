
  outName <- "new"
    
  library(magrittr)
  
  runFrom <- 0
  runTo <- 50
  excludes <- 34
  
  testing <- T
  testRmd <- F
  commitToGit <- T
  
  getNewData <- FALSE
  
  runReason <- paste0("Update to Rmd after overhaul of individual scripts. Testing = ",testing,".")
  
  dir() %>%
    grep("^\\d{4}_.*\\.R$",.,value=TRUE) %>%
    #grep(excludes,.,value = TRUE,invert = TRUE) %>%
    setNames(stringr::str_extract(.,"\\d{4}")) %>%
    `[` (names(.)[as.numeric(names(.)) <= runTo & as.numeric(names(.)) >= runFrom & !as.numeric(names(.)) %in% excludes]) %>%
    purrr::walk(source, verbose = TRUE)
  