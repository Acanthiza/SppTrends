
  outName <- "new"
    
  if(outName == "new") rm(list=grep("run",ls(),value=TRUE,invert = TRUE))
  
  
  runReason <- "Added full reporting rate and list length analyses. Testing that workflow"
  

  library(magrittr)
  
  runFrom <- 0
  runTo <- 60
  
  testing <- TRUE
  testRmd <- FALSE
  commitToGit <- TRUE
  
  getNewData <- FALSE
  
  dir() %>%
    grep("^\\d{4}_.*\\.R$",.,value=TRUE) %>%
    setNames(stringr::str_extract(.,"\\d{4}")) %>%
    `[` (names(.)[as.numeric(names(.)) <= runTo & as.numeric(names(.)) >= runFrom]) %>%
    purrr::walk(source, verbose = TRUE)
  