
  taxaModsOverall <- taxaModsOverall %>%
    dplyr::mutate(overall = pmap(list(Taxa
                                      , Common
                                      , resRR
                                      , resLL
                                      , resOcc
                                      )
                                 , function(a,b,c,d,e) year_difference_overall(a
                                                                             , b
                                                                             , c$yearDifferenceDf
                                                                             , d$yearDifferenceDf
                                                                             , e$yearDifferenceDf
                                                                             )
                                 )
                  )
  