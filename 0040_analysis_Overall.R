
  taxaMods <- ls(pattern = "taxaMods[[:upper:]]{2}") %>%
    tibble::enframe(name = NULL, value = "obj") %>%
    dplyr::mutate(data = map(obj,get)) %>%
    tidyr::unnest(cols = c(data)) %>%
    dplyr::mutate(yearDiffDfPath = map_chr(modPath,~gsub("Mod","YearDif",.))
                  , yearDiffDfPath = map_chr(yearDiffDfPath,~gsub("rds","feather",.))
                  )
  
  taxaModsOverall <- taxaMods %>%
    dplyr::mutate(yearDiffDfDo = map_lgl(yearDiffDfPath,file.exists)) %>%
    dplyr::filter(yearDiffDfDo) %>%
    dplyr::mutate(yearDiffDf  = map(yearDiffDfPath,read_feather)) %>%
    dplyr::select(type,yearDiffDf) %>%
    tidyr::unnest(cols = c(yearDiffDf)) %>%
    tidyr::nest(data = -c(Taxa,Common)) %>%
    dplyr::mutate(overall = pmap(list(Taxa
                                      , Common
                                      , data
                                      )
                                   , function(a,b,c) year_difference_overall(a
                                                                             , b
                                                                             , c
                                                                             )
                                 )
                  )
    