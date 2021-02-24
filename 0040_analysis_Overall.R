
  taxaModsOverall <- taxaModsOverall %>%
    dplyr::mutate(overall = pmap(list(Taxa
                                      , Common
                                      , resRR
                                      , resLL
                                      )
                                 , function(a,b,c,d) year_difference_overall(a
                                                                             , b
                                                                             , c$yearDifferenceDf
                                                                             , d$yearDifferenceDf
                                                                             )
                                 )
                  )
  