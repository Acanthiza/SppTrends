
  # connect to datamart
  con <- dbConnect(odbc::odbc()
                   , "BDBSA-via-DataMart"
                   , database = "BDBSA_DataMart"
                   , uid = read_csv("E:/key.csv") %>%
                     dplyr::filter(Secret == "BDBSA_user") %>%
                     pull(String)
                   , pwd = read_csv("E:/key.csv") %>%
                     dplyr::filter(Secret == "BDBSA_pwd") %>%
                     pull(String)
                   )
  
  # Link to each table
  # Survey
  sur <- dplyr::tbl(con, in_schema("Stage", "LUSURVEYNAME")) %>%
    dplyr::select(SURVEYNR,SURVEYNAME,SURVEYNAMEFULL)
  
  # Patch
  pat <- dplyr::tbl(con, in_schema("Stage", "SUPATCH")) %>%
    dplyr::select(SURVEYNR,PATCHID,LATITUDE,LONGITUDE,RELIABNR,LOCCOMM,EASTING,NORTHING) 
  
  # Visit
  vis <- dplyr::tbl(con, in_schema("Stage", "SUVISIT")) %>%
    dplyr::select(PATCHID,VISITNR,VISITDATE)
  
  # Species
  spp <- tbl(con, in_schema("Stage", "SUSPECIES")) %>%
    dplyr::filter(SPECIESTYPE != "P"
                  , DATEACCURACY != "C"
                  , DATEACCURACY != "T"
                  , ISCERTAIN == "Y"
                  , !NUMOBSERVED %in% c("0","none detected","None detected")
                  ) %>%
    dplyr::select(SPECIESTYPE,VISITNR,SPSEQNR,NSXCODE,SPECIESNR,OBSDATE,TIME,NUMOBSERVED,METHODNR,ISCERTAIN,SEX,DATEACCURACY,RECSTATUSCODE)
  
  luMethod <- tbl(con, in_schema("Stage","LUCOLLMETHOD")) %>%
    dplyr::select(1,2)
  
  # Taxaa lookup from 'not synonymous and not renamed' linked to to FL_FLSP
  luTaxa <- tbl(con, in_schema("Stage", "VSVNONSYNNOTREN"))
  
  # Get all patches
  patchesBDBSAAll <- sur %>%
    dplyr::right_join(pat) %>%
    dplyr::left_join(luRel, copy = TRUE) %>%
    dplyr::filter(maxDist <= setDist
                  , !is.na(LATITUDE)
                  ) %>%
    dplyr::collect()
  
  # Get all records
  taxaBDBSAAll <- spp %>%
    dplyr::left_join(vis, by = "VISITNR") %>%
    dplyr::left_join(luTaxa, by = "NSXCODE") %>%
    dplyr::left_join(luMethod) %>%
    dplyr::select(PATCHID
                  , VISITNR
                  , VISITDATE 
                  , SPECIES
                  , COMNAME1
                  , NSXCODE
                  , SPECIESNR = SPECIESNR.x
                  , ISINDIGENOUSFLAG
                  , METHODDESC 
                  ) %>%
    dplyr::collect() %>%
    dplyr::inner_join(patchesBDBSAAll, by = "PATCHID") %>%
    dplyr::filter(!grepl("in-active|sfossil",METHODDESC)
                  , !is.na(NSXCODE)
                  ) %>%
    dplyr::count(PATCHID,VISITDATE,SPECIES,NSXCODE,ISINDIGENOUSFLAG,SURVEYNR,SURVEYNAME,SURVEYNAMEFULL,LATITUDE,LONGITUDE,maxDist)
  
  write_rds(patchesBDBSAAll %>% dplyr::inner_join(taxaBDBSAAll %>% dplyr::count(PATCHID)), paste0(saveTo,"/patchesBDBSAAll.rds"))
  write_rds(taxaBDBSAAll, paste0(saveTo,"/taxaBDBSAAll.rds"))
  
  dbDisconnect(con)
