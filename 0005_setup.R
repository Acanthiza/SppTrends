
  library("timeR")
  timer <- createTimer()
  timer$start("setup")

#------Options-------
  
  # Reliability (distance) for accepting a patch
  setDist <- 1000
  
  # Buffer around polyMask in METRES
  polyBuffer <- 0
  
  # What is the area of interest (AOI) for this analysis?
  aoiName <- c("FLB","KAN")
  aoiFullName <- "MLR"
  Statewide <- F
  
  # Taxonomic grouping
  taxGroup <- "Order"
  
  # Length of one side of grid in metres
  gridSize <- 5000
  
  minRRListLength <- 6
  minLLListLength <- 2
  minListsPerYear <- # now set as a fraction of length(unique(cell or ibrasub or whatever))
  minTaxaOccurence <- 10
  
  minYear <- lubridate::year(Sys.time()) - 15
  
  
  quantProbs <- c(0.05, 0.5, 0.95)
  
  #---------Tidy-------
  # earliest and latest 'date' (= yearmon or YearMonth, i.e. 202012§)
  nowYearMon <- as.numeric(paste0(lubridate::year(Sys.Date()),sprintf("%02d",lubridate::month(Sys.Date()))))
  firstYearMon <- 183601
  
  # Which fields to collect from each data source - i.e. map the original source fields to these...
  collectFields <- c("LATITUDE","LONGITUDE"
                     ,"SPECIES", "CommonName"
                     ,"METHODDESC","NUMOBSERVED","ISINDIGENOUSFLAG","Rank"
                     ,"year","month"
                     ,"maxDist"
                     )
  
  nonRecords <- c("0","none detected","None detected","none detected ")

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
                       ,"feather"
                       ,"rgbif"
                       ,"knitr"
                       ,"bookdown"
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
                       ,"classInt"
                       ,"grid"
                       ,"gridExtra"
                       ,"maptools"
                       ,"cluster"
                       ,"alphahull"
                       ,"rangeBuilder"
                       ,"rstan"
                       ,"rstanarm"
                       ,"ggridges"
                       ,"parallelDist"
                       )
                     )
  
  purrr::walk(packages,library,character.only=TRUE)  
  
  
  
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
  
  
  #-------Load functions-------
  
  commonFiles <- path("..","template","toCommon")
  
  if(file.exists(commonFiles)){
    
    files <- dir_ls(commonFiles)
    newFiles <- files %>% gsub(commonFiles,path("common"),.)
    dir_create("common")
    walk2(files,newFiles,file_copy,overwrite=TRUE)
    
  }
  
  source("common/functions.R") # these are generic functions (e.g. vec_to_sentence)
  
  
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
  project <- basename(getwd())
  
  # What folder to save outputs to
  if(!exists("outName")) outName <- paste0(format(Sys.time(),"%Y-%m-%d-%H%M"),"_",aoiFullName)
  
  # Directory names
  outDir <- fs::path("out",outName)
  dir_create(outDir)
  

#-------Lookups------
  
  # LSAs
  luPolys <- tribble(
    ~LSA, ~REGION, ~LSARegion, ~Zone, ~System, ~R, ~G, ~B, ~A
    ,"HF","Hills and Fleurieu","Hills and Fleurieu",	"Agricultural",	"High-rainfall", 89, 23, 138, 255
    ,"AW","Alinytjara Wilurara","Alinytjara Wilurara",	"Arid",	"Arid",	191,	54,	44,	255
    ,"EP","Eyre Peninsula","Eyre Peninsula",	"Agricultural",	"Low-rainfall",	86,	156,	190,	255
    ,"KI","Kangaroo Island","Kangaroo Island",	"Agricultural",	"High-rainfall",	0,	137,	152,	255
    ,"NY","Northern and Yorke","Northern and Yorke",	"Agricultural",	"Low-rainfall",	248,	185,	44,	255
    ,"SAAL","South Australian Arid Lands","South Australian Arid Lands",	"Arid",	"Arid",	213,	94,	0,	255
    ,"MR","Murraylands and Riverland","Murraylands and Riverland",	"Agricultural",	"Low-rainfall",	0,	132,	197,	255
    ,"LC","Limestone Coast","Limestone Coast","Agricultural",	"High-rainfall",	102,	181,	98,	255
    ,"GA","Green Adelaide","Green Adelaide","Urban","High-rainfall",	34,	139,34,	255
    ) %>%
    dplyr::mutate(across(where(is.character),~gsub("Wilur|Wilu\\?",paste0("Wilu","\u1E5F"),.)))
  
  luRel <- read_csv("data/luReliability.csv")
  
  luRank <- tribble(
    ~Rank, ~sort
    , "Kingdom", 1
    , "Phylum", 2
    , "Class", 3
    , "Order", 4
    , "Family", 5
    , "Genus", 6
    , "Species", 7
    , "Subspecies", 8
    , "Variety", 9
    , "Form", 10
  )
  
  luLikelihood <- tribble(
    ~likelihood, ~maxVal
    , "Exceptionally unlikely", 0.01
    , "Extremely unlikely", 0.05
    , "Very unlikely", 0.1
    , "Unlikely", 1/3
    , "About as likely as not", 2/3
    , "Likely", 0.9
    , "Very likely", 0.95
    , "Extremely likely", 0.99
    , "Virtually certain", 1
    ) %>%
    dplyr::mutate(likelihood = fct_inorder(likelihood)
                  , range = cut(maxVal
                                   , breaks = c(0,.$maxVal)
                                   )
                  )

#----------Maps-----------
  
  # State map
  sa <- st_read("shp/aust_cd66states.shp"
                , crs = 4326
                , quiet = TRUE
                ) %>%
    dplyr::filter(STE == 4) %>%
    st_transform(crs = 3577)
  
  # Polygons
  # polys <- st_read("shp/LSA.shp"
  #                  , quiet = TRUE
  #                  ) %>%
  #   as_tibble() %>%
  #   dplyr::mutate(across(where(is.character),~gsub("Wilur|Wilu\\?",paste0("Wilu","\u1E5F"),.))) %>%
  #   dplyr::left_join(luPolys) %>%
  #   st_as_sf() %>%
  #   st_transform(crs = 3577)
  # 
  # polyPalette <- luPolys %>%
  #   dplyr::mutate(colour = rgb(R,G,B,A,maxColorValue = 255)) %>%
  #   dplyr::pull(colour, name = "REGION")
  
  # IBRA Sub
  ibraSub <- st_read("shp/LANDSCAPE_IbraSubregionAust70.shp") %>%
    st_transform(crs = 3577)

  polys <- ibraSub

#---------Parallel-----------
  
  # Cores to use for any parallel processing
  useCores <- if(parallel::detectCores() > 20) 20 else parallel::detectCores()-1
  
  # Plan for any furrr functions
  plan(multiprocess
       , workers = useCores
       , gc = TRUE
       )
  
  # For randomForestSRC
  options(rf.cores = useCores)
  
  
#---------References-------
  
  packageBibFile <- "packageRefs.bib"
  
  write_bib(packages
            ,file=packageBibFile
            ,tweak=TRUE
            ,width=1000
            )
  
  refs <- fix_bib(packageBibFile,isPackageBib = TRUE)

  
#-------
  
  timer$stop("setup")
  
  