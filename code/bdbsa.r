
  # Get all records
  outFile <- path("out","bdbsa","rawBDBSA.feather")
  
  get_all_bdbsa <- function(saveFile = outFile) {
    
    (tictoc::tic())
    
    # connect to BDBSA
    con <- dbConnect(odbc::odbc()
                     , "BDBSA Production"
                     , database = "BDBSA Production’"
                     , uid = Sys.getenv("BDBSA_PRD_user")
                     , pwd = Sys.getenv("BDBSA_PRD_pwd")
                     )
    
    # Link to each table
    excludeVars <- c("CREATED_DATE", "CREATED_USER", "MODIFIED_DATE", "MODIFIED_USER"
                     , "RECSTATUSCODE" # this is in patch and spp
                     )
    
    # Survey
    sur <- dplyr::tbl(con,"LUSURVEYNAME") %>%
      dplyr::select(!any_of(excludeVars))
    
    # Patch
    pat <- dplyr::tbl(con,"SUPATCH") %>%
      dplyr::select(!any_of(excludeVars)) %>%
      dplyr::filter(!is.na(LATITUDE))
    
    # Visit
    vis <- dplyr::tbl(con,"SUVISIT") %>%
      dplyr::select(!any_of(excludeVars))
    
    # Species
    spp <- tbl(con,"SUSPECIES") %>%
      dplyr::filter(SPECIESTYPE != "P"
                    , DATEACCURACY != "C"
                    , DATEACCURACY != "T"
                    , ISCERTAIN == "Y"
                    , !NUMOBSERVED %in% nonRecords
                    ) %>%
      dplyr::select(!any_of(excludeVars)) %>%
      dplyr::filter(!is.na(NSXCODE))
    
    # Method lookup
    luMethod <- tbl(con, "LUCOLLMETHOD")
    
    # Taxa lookup from 'not synonymous and not renamed' linked to to FL_FLSP
    luTaxa <- tbl(con, "VSVNONSYNNOTREN")
    
    rawBDBSA <- sur %>%
      dplyr::left_join(pat) %>%
      dplyr::left_join(vis) %>%
      dplyr::left_join(spp) %>%
      dplyr::left_join(luMethod) %>%
      dplyr::left_join(luTaxa, by = "NSXCODE") %>%
      dplyr::collect()
    
    write_feather(rawBDBSA,saveFile)
    
    dbDisconnect(con)
    
    (tictoc::toc())
    
    return(rawBDBSA)
    
  }
  
  rawBDBSA <- if(getNewData|!file.exists(outFile)) get_all_bdbsa() else read_feather(outFile)
