

  library("bookdown")
  library("tidyverse")
  library("fs")
  
  if(!exists("outName")) outName <- "2021-01-13-1158_KI"
  
  # Ensure any previous failed knit result is cleaned up
  file_delete("_main.Rmd")
  
  # Make book, collecting all .Rmd files in order and knitting to book
  render_book("Rmd")
  
  #rmarkdown::render("_test.Rmd")
  
  
  