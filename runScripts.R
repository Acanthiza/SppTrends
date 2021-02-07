
  outName <- "2021-02-04-0653_MLR"
    
  if(outName == "new") rm(list=grep("run",ls(),value=TRUE,invert = TRUE))
  
  
  runReason <- "First attempt to run all taxa through rr and ll."
  

  library(magrittr)
  
  runFrom <- 0
  runTo <- 30
  
  testing <- F
  testRmd <- FALSE
  commitToGit <- F
  
  getNewData <- FALSE
  
  dir() %>%
    grep("^\\d{4}_.*\\.R$",.,value=TRUE) %>%
    setNames(stringr::str_extract(.,"\\d{4}")) %>%
    `[` (names(.)[as.numeric(names(.)) <= runTo & as.numeric(names(.)) >= runFrom]) %>%
    purrr::walk(source, verbose = TRUE)
  