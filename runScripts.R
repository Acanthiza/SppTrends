
  outName <- "new"
    
  library(magrittr)
  
  runFrom <- 30
  runTo <- 50
  excludes <- 35
  
  testing <- T
  testRmd <- F
  commitToGit <- F
  
  getNewData <- FALSE
  
  runReason <- paste0("'Change' is now comparison between a reference year and a recent year. Testing = ",testing,".")
  
  dir() %>%
    grep("^\\d{4}_.*\\.R$",.,value=TRUE) %>%
    #grep(excludes,.,value = TRUE,invert = TRUE) %>%
    setNames(stringr::str_extract(.,"\\d{4}")) %>%
    `[` (names(.)[as.numeric(names(.)) <= runTo & as.numeric(names(.)) >= runFrom & !as.numeric(names(.)) %in% excludes]) %>%
    purrr::walk(source, verbose = TRUE)
  