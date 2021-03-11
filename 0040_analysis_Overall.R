
  taxaMods <- ls(pattern = "taxaMods[[:upper:]]{2}") %>%
    tibble::enframe(name = NULL, value = "obj") %>%
    dplyr::mutate(data = map(obj,get)) %>%
    tidyr::unnest(cols = c(data))
  
  taxaModsOverall <- taxaMods %>%
    dplyr::mutate(yearDiffDf = map(res,"yearDifferenceDf")) %>%
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
    