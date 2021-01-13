#------Options-------

  # CLEAN RUN - set to TRUE for rebuilding all objects
  cleanRun <- F
  
  # Run which chunks...
  bdbsa <- T
  
  trend <- T
  
  # Realiability (distance) for accepting a patch
  setDist <- 1000
  
  # Buffer around polyMask
  polyBuffer <- 15000
  
  # Mask used to select PATCHIDs
  polyMask <- c(NULL
                , "Adelaide & Mt Lofty Ranges"
                , "South Australian Arid Lands"
                , "Alinytjara Wilurara"
                , "South Australian Murray-Darling Basin"
                , "Northern & Yorke"
                , "Eyre Peninsula"
                , "Kangaroo Island"
                , "South East"
                )
  
  # maximium number of levels to plot
  maxLevels <- 25
  
  # Length of one side of grid in metres
  gridSize <- 10000
  
  minRRListLength <- 6
  minLLListLength <- 2
  minListsPerYear <- # now set as a fraction of length(unique(cell or ibrasub or whatever))
  minTaxaOccurence <- 4
  
  # Stan chains options
  # quick
  # doChains <- 4
  # doIter <- 1000
  
  # good
  doChains <- 5
  doIter <- 3000

#-----Packages-----
  
  packages <- unique(c("base"
                       ,"knitr"
                       ,"bookdown"
                       ,"readr"
                       ,"dplyr"
                       ,"tidyr"
                       ,"purrr"
                       ,"furrr"
                       ,"tibble"
                       ,"stringr"
                       ,"forcats"
                       ,"ggplot2"
                       ,"GGally"
                       ,"lubridate"
                       ,"dbplyr"
                       ,"DBI"
                       ,"fs"
                       ,"sf"
                       ,"raster"
                       ,"snowfall"
                       ,"corrplot"
                       ,"caret"
                       ,"rasterVis"
                       ,"tabularaster"
                       ,"parallel"
                       ,"doParallel"
                       ,"ggrepel"
                       ,"gridExtra"
                       ,"DT"
                       ,"tmap"
                       ,"DataExplorer"
                       ,"mgcv" # referenced in geom_smooth
                       ,"classInt"
                       ,"grid"
                       ,"gridExtra"
                       ,"maptools"
                       ,"cluster"
                       ,"alphahull"
                       ,"rangeBuilder"
                       ,"kableExtra"
                       ,"rstan"
                       ,"rstanarm"
                       ,"ggridges"
                       ,"parallelDist"
                       )
                     )
  
  purrr::walk(packages,library,character.only=TRUE)  
  
  write_bib(packages,file="common/packageCitations.bib",tweak=TRUE)
  
#------Options-------
  
  # rstan options
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
  
  
  # tmap options
  tmap::tmap_options(basemaps = c("OpenStreetMap.Mapnik"
                                  , "Esri.WorldImagery"
                                  )
                     )
  
  tmap::tmap_mode("view")
  
  
  # Scientific notation?
  options(scipen = 999)
  
  
  # ggplot theme
  ggplot2::theme_set(ggplot2::theme_grey())
  
  
#-------Get functions-------
  
  if(file.exists("../Tools/functions.R")) file_copy("../Tools/functions.R","common/functions.R",overwrite = TRUE)
  source("common/functions.R")
  
  if(file.exists("../Tools/getGBIFTax.R")) file_copy("../Tools/getGBIFTax.R","common/getGBIFTax.R",overwrite = TRUE)
  source("common/getGBIFTax.R")
  
  if(file.exists("../../git/Refs/refs.bib")) file_copy("../../git/Refs/refs.bib","common/refs.bib",overwrite = TRUE)
  
#------Chunk options-------
  # Set default chunk options (can adjust individual chunks differently if required)
  
  knitr::opts_chunk$set(echo = FALSE
                        , warning = FALSE
                        , error = FALSE
                        , message = FALSE
                        , comment = NA
                        )
  
  options(knitr.kable.NA = ""
          , kableExtra.html.bsTable = TRUE
          )
  
  
#------Project------
  
  # Project name
  project <- "SppTrends"
  
  # What folder (= 'run') to save outputs to
  lastRun <- dir_info(paste0("../../../Data/",project)) %>%
    dplyr::filter(type == "directory") %>%
    dplyr::select(path) %>%
    dplyr::mutate(run = as.numeric(map_chr(path,str_extract,"[[:digit:]]+"))) %>%
    pull(run) %>%
    max(na.rm=TRUE)
  
  cleanRun <- if(exists("cleanRun")) cleanRun else FALSE
  
  run <- if(cleanRun) {
    
    scales::ordinal(lastRun + 1)
    
  } else {
    
    scales::ordinal(lastRun)
    
  }
  
  # Directory names
  outDir <- paste0("../../../Data/",project,"/",run)
  rasterFolder <- "C:/Users/nwilloughby/Documents/Data/Raster/Env_Cover_Layers_30m"
  
  # Make directories
  dir_create(outDir)

#-------Lookups------
  
  luRel <- read_csv("data/luReliability.csv")
  
  

#----------Maps-----------
  
  # State map
  sa <- st_read("shp/aust_cd66states.shp", crs = 4326) %>%
    dplyr::filter(STE == 4)
  
  # Polygons
  polys <- st_read("shp/nrm.shp")
  
  # IBRA Sub
  ibraSub <- st_read("shp/LANDSCAPE_IbraSubregionAust70.shp")


#---------Parallel-----------
  
  # Cores to use for any parallel processing
  useCores <- parallel::detectCores()/2
  
  # Plan for any furrr functions
  plan(multiprocess
       , workers = useCores
       , gc = TRUE
       )
  
  # For randomForestSRC
  options(rf.cores = useCores)