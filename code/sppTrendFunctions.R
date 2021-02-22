  
  filter_taxa_data <- function(df
                               , minListLengthThresh = 3
                               , maxListLengthOccurenceThresh = 3
                               , minlistOccurenceThresh = 5
                               , minYearsThresh = 3
                               ) {
    
    find_min_list_length <- function(df) {
      
      min(df$listLength)
      
    }
    
    find_max_list_length_occurence <- function(df) {
      
      df %>%
        dplyr::count(geo2,Taxa,listLength,name="lists") %>%
        dplyr::group_by(Taxa) %>%
        dplyr::filter(listLength == max(listLength)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(listLength == min(listLength)) %>%
        dplyr::pull(listLength) %>%
        unique()
      
    }
    
    find_min_list_occurence <- function(df) {
      
      df %>%
        dplyr::count(geo2,!!ensym(taxGroup),Taxa,name = "lists") %>%
        dplyr::filter(lists == min(lists)) %>%
        dplyr::pull(lists) %>%
        unique()
      
    }
    
    find_min_years <- function(df) {
      
      df %>%
        dplyr::count(geo2,!!ensym(taxGroup),Taxa,year,name = "blah") %>%
        dplyr::count(geo2,!!ensym(taxGroup),Taxa, name = "years") %>%
        dplyr::filter(years == min(years)) %>%
        dplyr::pull(years) %>%
        unique()
      
    }
    
    minListLength <- find_min_list_length(df)
    maxListLengthOccurence <- find_max_list_length_occurence(df)
    minlistOccurence <- find_min_list_occurence(df)
    minYears <- find_min_years(df)
    
    
    while(minListLength < minListLengthThresh |
          maxListLengthOccurence < maxListLengthOccurenceThresh |
          minlistOccurence < minlistOccurenceThresh |
          minYears < minYearsThresh
    ) {
      
      df <- df %>%
        dplyr::filter(listLength > minListLengthThresh)
      
      removeTaxaOnShortLists <- df %>%
        dplyr::count(geo2,!!ensym(taxGroup),Taxa,listLength,name="lists") %>%
        dplyr::group_by(Taxa) %>%
        dplyr::filter(listLength == max(listLength)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(listLength < maxListLengthOccurenceThresh) %>%
        dplyr::distinct(geo2,!!ensym(taxGroup),Taxa)
      
      removeTaxaWithFewOccurrences <- df %>%
        dplyr::anti_join(removeTaxaOnShortLists) %>%
        dplyr::count(geo2,!!ensym(taxGroup),Taxa,name = "lists") %>%
        dplyr::filter(lists < minlistOccurenceThresh) %>%
        dplyr::distinct(geo2,!!ensym(taxGroup),Taxa)
      
      removeTaxaWithFewYears <- df %>%
        dplyr::anti_join(removeTaxaOnShortLists) %>%
        dplyr::anti_join(removeTaxaWithFewOccurrences) %>%
        dplyr::count(geo2,!!ensym(taxGroup),Taxa,year,name = "blah") %>%
        dplyr::count(geo2,!!ensym(taxGroup),Taxa, name = "years") %>%
        dplyr::filter(years < minYearsThresh) %>%
        dplyr::distinct(geo2,!!ensym(taxGroup),Taxa)
      
      df <- df %>%
        dplyr::anti_join(removeTaxaOnShortLists) %>%
        dplyr::anti_join(removeTaxaWithFewOccurrences) %>%
        dplyr::anti_join(removeTaxaWithFewYears) %>%
        dplyr::add_count(list, name = "listLength")
      
      minListLength <- find_min_list_length(df)
      maxListLengthOccurence <- find_max_list_length_occurence(df)
      minlistOccurence <- find_min_list_occurence(df)
      minYears <- find_min_years(df)
      
    }
    
    return(df)
    
  }
  