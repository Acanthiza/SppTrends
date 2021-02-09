
  outName <- "2021-02-07-1313_MLR"
    
  if(outName == "new") rm(list=grep("run",ls(),value=TRUE,invert = TRUE))
  
  
  runReason <- "Fixed issue with wrong object in `taxaSummary.rmd`."
  

  library(magrittr)
  
  runFrom <- 0
  runTo <- 50
  
  testing <- T
  testRmd <- F
  commitToGit <- T
  
  getNewData <- FALSE
  
  dir() %>%
    grep("^\\d{4}_.*\\.R$",.,value=TRUE) %>%
    setNames(stringr::str_extract(.,"\\d{4}")) %>%
    `[` (names(.)[as.numeric(names(.)) <= runTo & as.numeric(names(.)) >= runFrom]) %>%
    purrr::walk(source, verbose = TRUE)
  