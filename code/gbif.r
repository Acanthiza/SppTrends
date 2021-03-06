
  # Get all records
  outFile <- path("out","gbif","rawGBIF.feather")
  fs::dir_create(dirname(outFile))
  
  #--------Get gbif data----------
  if(getNewData|!file.exists(outFile)){
  
    aoiWKT <- st_as_sfc(st_bbox(sa)) %>%
      st_buffer(100000) %>%
      st_geometry() %>%
      st_transform(crs = 4326) %>%
      st_as_text()
    
    temp <- occ_download(pred("taxonKey",1) # 1 = Animalia
                      , pred("hasCoordinate",TRUE)
                      , pred("geometry",aoiWKT)
                      , user = Sys.getenv("GBIF_user")
                      , pwd = Sys.getenv("GBIF_pwd")
                      , email = Sys.getenv("GBIF_email")
                      )
    
    occ_download_wait(temp)
    
    meta <- occ_download_meta(temp)
    
    getDownload <- occ_download_get(temp
                                    , path = path("out","gbif")
                                    , overwrite = TRUE
                                    )
    
  }
  
  
  #------get gbif meta data---------
  
  metaGBIF <- list()
    
  info <- occ_download_list(Sys.getenv("GBIF_user"),Sys.getenv("GBIF_pwd"))$results %>%
    dplyr::mutate(created = ymd_hms(created,tz=Sys.timezone())) %>%
    dplyr::filter(created == max(created))
  
  metaGBIF$key <- info %>%
    dplyr::pull(key)
  
  metaGBIF$doi <- info %>%
    dplyr::pull(doi)
  
  metaGBIF$licence <- info %>%
    dplyr::pull(license)
  
  metaGBIF$date <- info %>%
    dplyr::pull(created)
  
  metaGBIF$ref <- get_bib(metaGBIF$doi,outFile = "gbifDataRef.bib")
  
  # Make a reference for the download
  ref <- read_lines("gbifDataRef.bib")
  ref[1] <- paste0("@misc{GBIFRef,")
  write_lines(ref,"gbifDataRef.bib")
  
  
  #-------unzip gbif data--------
  
  outFile <- path("out","gbif",metaGBIF$key,"occurrence.txt")
  
  if(!file.exists(outFile)) {
    
    unzip(path("out","gbif",paste0(metaGBIF$key,".zip"))
          , exdir = path("out","gbif",metaGBIF$key)
          )
    
    rawGBIF <- data.table::fread(outFile)
    
    write_feather(rawGBIF,path("out","gbif","rawGBIF.feather"))
    
    # Get original citations for downloaded data
    
    # get_dataset_DOI <- function(xmlFile) {
    #   
    #   XML::xmlParse(xmlFile) %>%
    #     XML::getNodeSet("//*/dataset") %>%
    #     XML::xmlToDataFrame() %>%
    #     dplyr::pull(alternateIdentifier)
    #   
    # }
    
    # datasetRefs <- dir_info(path("out","gbif",metaGBIF$key,"dataset")) %>%
    #   dplyr::select(path) %>%
    #   dplyr::mutate(DOI = map_chr(path,get_dataset_DOI))
    # 
    # unlink("datasetRefs.bib")
    # 
    # future_walk(c(metaGBIF$doi,datasetRefs$DOI)
    #             ,get_bib
    #             ,outFile = "datasetRefs.bib"
    #             )
    # 
    # fix_bib("datasetRefs.bib"
    #         ,isPackageBib = FALSE
    #         ,makeKey = TRUE
    #         )
    
  }
  
  
  #------load data-------
  
  rawGBIF <- read_feather(path("out","gbif","rawGBIF.feather"))
  
  # refsGBIF <- bib2df::bib2df("datasetRefs.bib")
  
  
  