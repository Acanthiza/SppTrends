
  
#-------make book-------
  
  write_feather(getTimer(timer),path(outDir,"timer.feather"))
  
  source("makeBook.R")
  
#-------Export---------
  
  toPath <- path("//env.sa.gov.au/dfsroot/IST/DEHProjects/Landscapes","SppTrends",outName)
  dir_create(toPath)
  
  if(file.exists(toPath)){
    
    file.copy(path("_book","_main.docx")
              , path(toPath,"Report.docx")
              , overwrite = TRUE
              )
    
  }
  