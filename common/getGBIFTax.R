
  
# required packages

  library("rgbif")
  
  
# a function to retrieve taxonomy to accepted names and retrieve taxonomic hierarchy for a df with a column of taxonomic names

  gbif_tax <- function(df,sppCol=1,idCol=2,outFile="data/luGBIF.csv",kingType="Plantae"){
    
    sppCol <- if(is.character(sppCol)) sppCol else names(df)[sppCol]
    idCol <- if(is.character(idCol)) idCol else names(df)[idCol]
    
    dat <- if(file.exists(outFile)) {
      df %>%
        dplyr::rename(Species = sppCol
                      , id = idCol
                      ) %>%
        dplyr::anti_join(read_csv(outFile) %>%
                           dplyr::filter(!is.na(key))
                        , by = "id"
                         ) %>%
        dplyr::arrange(Species)
      } else if (file.exists(paste0(gsub(".csv","",outFile),"_temp.csv"))) {
        df %>%
          dplyr::rename(Species = sppCol
                        , id = idCol
                        ) %>%
          dplyr::anti_join(read_csv(paste0(gsub(".csv","",outFile),"_temp.csv"))) %>%
          dplyr::arrange(Species)
        } else {
          df %>%
            dplyr::rename(Species = sppCol
                          , id = idCol
                          ) %>%
            dplyr::arrange(Species)
          }
    
    taxa <- dat %>%
      pull(Species) %>%
      gsub("\\s*\\(.*\\)|\\'|dead|\\?| spp\\.| sp\\.|#","",.) %>%
      gsub(" x .*$| X .*$","",.)
    
    resStructure <- tribble(~id
                            , ~usageKey
                            , ~scientificName
                            , ~canonicalName
                            , ~rank
                            , ~status
                            , ~confidence
                            , ~matchType
                            , ~kingdom
                            , ~phylum
                            , ~class
                            , ~order
                            , ~family
                            , ~genus
                            , ~species
                            , ~common
                            , ~kingdomKey
                            , ~phylumKey
                            , ~classKey
                            , ~orderKey
                            , ~familyKey
                            , ~genusKey
                            , ~speciesKey
                            , ~Species
                            , ~Stamp
                            ) %>%
        dplyr::mutate_all(as.character())
    
    
    if(length(taxa)>0){
      
      for (i in 1:length(taxa)){
        
        print(taxa[i])
        
        taxGBIF <- name_backbone(taxa[i], kingdom = kingType)
        
        
        if(taxGBIF$matchType == "NONE") {
          
          taxGBIF <- tibble(Species = taxa[i])
          
        }
        
        
        taxGBIF <- if(sum(grepl("acceptedUsageKey",names(taxGBIF)))>0) {
          name_usage(taxGBIF$acceptedUsageKey,return="data") %>%
            dplyr::mutate(matchType = "Synonym") %>%
            dplyr::rename(usageKey = key
                          , status = taxonomicStatus
                          )
          } else {
            taxGBIF
          }
        
        comGBIF <- name_usage(taxGBIF$usageKey, language = "ENGLISH", data="vernacularNames")$data %>%
            dplyr::filter(language == "eng")
        
        taxGBIF[["common"]] <- if(nrow(comGBIF) > 0) {
          
          comGBIF %>%
            dplyr::count(vernacularName) %>%
            dplyr::arrange(desc(n)) %>%
            dplyr::slice(1) %>%
            dplyr::pull(vernacularName)
          
          } else ""
        
        res <- resStructure %>%
          dplyr::bind_rows(
            if(exists("taxGBIF")) {
              as_tibble(taxGBIF) %>%
                dplyr::bind_cols(dat[i,]) %>%
                dplyr::mutate(Stamp = Sys.time()) %>%
                dplyr::mutate_all(as.character)
            } else {
              dat[i,] %>%
                dplyr::mutate(Stamp = Sys.time()) %>%
                dplyr::mutate_all(as.character)
            }
          ) %>%
          dplyr::select(1:ncol(resStructure))
          
        if(file.exists(paste0(gsub(".csv","",outFile),"_temp.csv"))) {
          
          write_csv(res
                    , paste0(gsub(".csv","",outFile),"_temp.csv")
                    , append = TRUE
                    )
          
        } else {
          
          write_csv(res
                    , paste0(gsub(".csv","",outFile),"_temp.csv")
                    )  
          
        }
        
        }
      
      }
    
    # Clean up results
    if(file.exists(paste0(gsub(".csv","",outFile),"_temp.csv"))) {
      
      read_csv(paste0(gsub(".csv","",outFile),"_temp.csv")) %>%
        dplyr::mutate_if(is.logical,as.character) %>%
        dplyr::mutate(Taxa = if_else(!is.na(species)
                                        ,species
                                        ,if_else(!is.na(genus)
                                                 , genus
                                                 , if_else(!is.na(family)
                                                           , family
                                                           , if_else(!is.na(order)
                                                                     , order
                                                                     , if_else(!is.na(class)
                                                                               , class
                                                                               , if_else(!is.na(phylum)
                                                                                         , phylum
                                                                                         , if_else(!is.na(kingdom)
                                                                                                   , kingdom
                                                                                                   , Species
                                                                                                   )
                                                                                         )
                                                                               )
                                                                     )
                                                           )
                                                 )
                                        )
                      , Stamp = as.POSIXct(Stamp)
                      ) %>%
        dplyr::select(id = id
                      , key = usageKey
                      , Rank = rank
                      , Taxa
                      , originalName = Species
                      , Kingdom = kingdom
                      , Phylum = phylum
                      , Class = class
                      , Order = order
                      , Family = family
                      , Genus = genus
                      , Species = species
                      , Common = common
                      , scientificName
                      , canonicalName
                      , Status = status
                      , Confidence = confidence
                      , Match = matchType
                      , Stamp
                      ) %>%
        (if(!file.exists(outFile)) {
          
          function(x) x
          
        } else {
          
          function(x) x %>%
            dplyr::bind_rows(read_csv(outFile) #%>% dplyr::mutate(key = as.character(key)
                                                                 #, Confidence = as.character(key)
                                                                # )
                             ) %>%
            dplyr::mutate(Stamp = as.POSIXct(Stamp)) %>%
            dplyr::mutate_if(is.numeric,as.character)
          
        }
        
        ) %>%
        dplyr::group_by(id) %>%
        dplyr::arrange(desc(Stamp)) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        write_csv(outFile)
      
      file.remove(paste0(gsub(".csv","",outFile),"_temp.csv"))
      
    } else {
      
      {warning( "No taxa supplied" ) }
  
    }

  }  
  
  