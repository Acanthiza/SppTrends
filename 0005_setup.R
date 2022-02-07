
  library("timeR")
  timer <- createTimer()
  timer$start("setup")
  
#------Options-------
  
  # Reliability (distance) for accepting a patch
  setDist <- 1000
  
  # Buffer around polyMask in METRES
  polyBuf <- 1000
  
  # What is the area of interest (AOI) for this analysis?
  aoiName <- "Flinders Lofty Block"
  aoiFullName <- "Flinders Lofty Block"
  Statewide <- F
  
  # Which polygons define the AOI?
  polyMask <- c(NULL
                , "Simpson Strzelecki Dunefields"
                , "Stony Plains"
                , "Naracoorte Coastal Plain"
                , "Nullarbor"
                , "Southern Volcanic Plain"
                , "Riverina"
                , "Central Ranges"
                , "Murray Darling Depression"
                , "Flinders Lofty Block"
                , "Hampton"
                , "Kanmantoo"
                , "Channel Country"
                , "Great Victoria Desert"
                , "Broken Hill Complex"
                , "Finke"
                , "Eyre Yorke Block"
                , "Gawler"
                )
  
  # polyMask <- c(NULL
  #               , "Kangaroo Island"
  #               )
  
  # Geo context
  geo1 <- "IBRA_REG_N"
  geo2 <- "IBRA_SUB_N"
  
  # Taxonomic grouping
  taxGroup <- "class"
  
  # Length of one side of grid in metres
  innerGrid <- 1000
  outerGrid <- 10*innerGrid
  
  minYear <- 1990 #lubridate::year(Sys.time()) - 30
  
  quantProbs <- c(0.05, 0.5, 0.95)
  
  #---------Tidy-------
  
  # Which fields to collect from each data source - i.e. map the original source fields to these...
  collectFields <- c("LATITUDE","LONGITUDE"
                     ,"date","year","month","yearmon","yday"
                     ,"SPECIES", "CommonName"
                     ,"METHODDESC","NUMOBSERVED","ISINDIGENOUSFLAG","Rank"                     
                     ,"maxDist"
                     )
  
  nonRecords <- c("0","none detected","None detected","none detected ")
  
  
  #------Settings - cluster------
  
  # From the fastcluster manual https://cran.r-project.org/web/packages/fastcluster/vignettes/fastcluster.pdf
  # "The following three methods are intended for Euclidean data only, ie. when X contains
  # the pairwise squared distances between vectors in Euclidean space. The algorithm
  # will work on any input, however, and it is up to the user to make sure that applying
  # the methods makes sense:
  # ....
  # * centroid
  # * median
  # * ward.D and Ward.D2"
  
  # Thus, those methods are not used...
  
  clustMethod <- tibble::tibble( method = c(
    "single" # had very few clusterings with minGroupSize > 1
    , "complete" # consistently well below par silhouette widths
    , "average" # didn't have any clusterings with minGroupSize > 1
    , "mcquitty"# had very few clusterings with minGroupSize > 1
    #, "ward.D"
    , "ward.D2" # Let this one back in despite the fastcluster quote above (and hclustgeo uses it)
    #, "centroid" # didn't have any clusterings with minGroupSize > 1
    #, "median" # had very few clusterings with minGroupSize > 1
  )
  )
  
  # Possible groups
  possibleGroups <- 2:20
  minSites <- 10
  minSpp <- 3
  
  #--------Analysis - RR----------------
  
  # Years at which to predict (and compare change)
  testYears <- tibble::tribble(~type, ~year
                       , "reference", 1995
                       , "recent", 2020
                       ) %>%
    tidyr::unnest(cols = c(year))
  
  reference <- testYears$year[testYears$type == "reference"]
  recent <- testYears$year[testYears$type == "recent"]

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
                       ,"gamm4"
                       ,"ggridges"
                       ,"parallelDist"
                       ,"ubms"
                       ,"red"
                       ,"cooccur"
                       )
                     )
  
  purrr::walk(packages,library,character.only=TRUE)  
  
  
  
  #-------Load functions-------
  
  commonFiles <- path("..","template","toCommon")
  
  if(file.exists(commonFiles)){
    
    files <- dir_ls(commonFiles)
    newFiles <- files %>% gsub(commonFiles,path("common"),.)
    dir_create("common")
    walk2(files,newFiles,file_copy,overwrite=TRUE)
    
  }
  
  source("common/functions.R") # these are generic functions (e.g. vec_to_sentence)
  source("code/sppTrendFunctions.R")
  
  
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
  outNameRaw <- outName
  if(outName == "new") outName <- paste0(format(Sys.time(),"%Y-%m-%d-%H%M"),"_",aoiFullName)
  
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
  
  
  # Diagnostics
  diagnostics <- tribble(
    ~diagnostic, ~diagDefinition
    , "avClustSize", "Average cluster size"
    , "maxClustSize", "Maximum cluster size"
    , "propClustersMin", paste0("Proportion of clusters with more than ",minSites," sites")
    , "propSitesClustersMin", paste0("Proportion of sites in clusters with more than ",minSites," sites")
    , "macroSil", "Mean silhouette width of clustering"
    , "macroWSS",  "Total within group sum-of-squares across all clusters"
    , "propGoodClusters", "Proportion of good clusters"
    , "propGoodSites", "Proportion of sites in good clusters"
    , "kappa", "Kappa value"
    , "accuracy", "Accuracy value"
  ) %>%
    dplyr::left_join(
      tribble(
        ~diagnostic, ~highGood
        , "avClustSize", TRUE
        , "maxClustSize", FALSE
        , "propClustersMin", TRUE
        , "propSitesClustersMin", TRUE
        , "macroSil", TRUE
        , "macroWSS",  FALSE
        , "propGoodClusters", TRUE
        , "propGoodSites", TRUE
        , "kappa", TRUE
        , "accuracy", TRUE
      )
    ) %>%
    dplyr::left_join(
      tribble(
        ~diagnostic, ~weightClusters
        , "avClustSize", FALSE
        , "maxClustSize", FALSE
        , "propClustersMin", TRUE
        , "propSitesClustersMin", TRUE
        , "macroSil", TRUE
        , "macroWSS",  TRUE
        , "propGoodClusters", FALSE
        , "propGoodSites", FALSE
        , "kappa", FALSE
        , "accuracy", FALSE
      )
    ) %>%
    dplyr::mutate(across(where(is.character),~fct_inorder(.)))

#----------Maps-----------
  
  # State map
  sa <- st_read("shp/aust_cd66states.shp"
                , crs = 4326
                , quiet = TRUE
                ) %>%
    dplyr::filter(STE == 4) %>%
    st_transform(crs = 3577)
  
  # Polygons
  LSA <- st_read("shp/LSA.shp"
                   , quiet = TRUE
                   ) %>%
    as_tibble() %>%
    dplyr::mutate(across(where(is.character),~gsub("Wilur|Wilu\\?",paste0("Wilu","\u1E5F"),.))) %>%
    dplyr::left_join(luPolys) %>%
    st_as_sf() %>%
    st_transform(crs = 3577)

  lsaPalette <- luPolys %>%
    dplyr::mutate(colour = rgb(R,G,B,A,maxColorValue = 255)) %>%
    dplyr::pull(colour, name = "REGION")
  
  # IBRA Sub
  ibraSub <- st_read("shp/LANDSCAPE_IbraSubregionAust70.shp") %>%
    st_transform(crs = 3577) %>%
    sf::st_make_valid() %>%
    st_intersection(sa)
  
  luGeo <- ibraSub %>%
    st_set_geometry(NULL) %>%
    dplyr::distinct(IBRA_SUB_N,IBRA_SUB_C,IBRA_SUB_1,IBRA_REG_N,IBRA_REG_C,IBRA_REG_1)
    

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
  
  
  
#------Options-------
  
  #------ rstan options-------
  options(mc.cores = useCores)
  rstan_options(auto_write = TRUE)
  Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
  
  testChains <- 3
  testIter <- 1000
  
  useChains <- 4
  useIter <- 2000
  
  
  # tmap options
  tmap::tmap_options(basemaps = c("OpenStreetMap.Mapnik"
                                  , "Esri.WorldImagery"
                                  )
                     )
  
  tmap::tmap_mode("plot")
  
  
  # Scientific notation?
  options(scipen = 999)
  
  
  # ggplot theme
  ggplot2::theme_set(ggplot2::theme_grey())
  
  
#---------References-------
  
  pac <- .packages()
  knitr::write_bib(pac, "package_citations.bib")
  refs <- bib2df::bib2df("package_citations.bib")
  
  
#------Settings - save or load -------
  
  settings <- fs::path(outDir,".RData")
  
  if(outNameRaw == "new") save.image(file = settings) #else load(settings)

  
#-------
  
  timer$stop("setup")
  
  