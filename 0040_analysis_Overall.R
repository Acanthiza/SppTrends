
  taxaModsOverall <- taxaMods %>%
    dplyr::full_join(taxaModsLL %>%
                       dplyr::select(Taxa,where(is.list))
                     , by = c("Taxa")
                     #, suffix = c("RR","LL")
                     ) %>%
    dplyr::full_join(taxaModsOcc %>%
                       dplyr::select(Taxa,where(is.list))
                     , by = c("Taxa")
                     )
  
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
    