

  library("bookdown")
  library("tidyverse")
  library("fs")
  
  if(!exists("outName")) outName <- "2021-01-13-1158_KI"
  
  # Ensure any previous failed knit result is cleaned up
  unlink("_main.Rmd")
  
  # Make book, collecting all .Rmd files in order and knitting to book
  render_book("Report.Rmd")
  
  #rmarkdown::render("_test.Rmd")
  
  
  