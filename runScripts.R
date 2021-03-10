
  outName <- "2021-03-10-0637_agricultural zone"
    
  library(magrittr)
  
  runFrom <- 0
  runTo <- 34
  excludes <- 0
  
  testing <- T
  testRmd <- F
  commitToGit <- F
  
  getNewData <- FALSE
  
  runReason <- paste0("Left original date (rather than lumping to month) in taxaAll. Testing = ",testing,".")
  
  dir() %>%
    grep("^\\d{4}_.*\\.R$",.,value=TRUE) %>%
    #grep(excludes,.,value = TRUE,invert = TRUE) %>%
    setNames(stringr::str_extract(.,"\\d{4}")) %>%
    `[` (names(.)[as.numeric(names(.)) <= runTo & as.numeric(names(.)) >= runFrom & !as.numeric(names(.)) %in% excludes]) %>%
    purrr::walk(source, verbose = TRUE)
  