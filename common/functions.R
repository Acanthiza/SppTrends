
#------Clustering---------

  make_clusters <- function(data
                            ,methodsDf = clustMethod
                            ,sppCol = "Taxa"
                            ,siteCol = "list"
                            ,groups = possibleGroups[possibleGroups < 0.5*length(unique(data$Taxa))]
                            ,minTaxaCount = 1
                            ) {
    
    datWide <- data %>%
      dplyr::add_count(!!ensym(sppCol)) %>%
      dplyr::filter(n > minTaxaCount) %>%
      dplyr::select(!!ensym(sppCol),!!ensym(siteCol)) %>%
      dplyr::mutate(p = 1) %>%
      tidyr::pivot_wider(names_from = all_of(sppCol), values_from = "p", values_fill = 0)
    
    siteNames <- datWide %>% dplyr::pull(!!ensym(siteCol))
    
    dist <- parDist(datWide %>% tibble::column_to_rownames(siteCol) %>% as.matrix()
                    , method = "bray"
                    , threads = useCores
                    )
    
    assign("sqDist",as.matrix(dist^2),pos = .GlobalEnv)
    
    dend <- methodsDf %>%
      dplyr::mutate(dend = map(method
                               ,~fastcluster::hclust(dist, .)
                               )
                    )
    
    clust <- dend %>%
      dplyr::mutate(clusters = map(dend
                                   , cutree
                                   , groups
                                   )
                    , clusters = map(clusters
                                     , as_tibble
                                     )
                    ) %>%
      dplyr::select(-dend) %>%
      tidyr::unnest(clusters) %>%
      tidyr::pivot_longer(2:ncol(.),names_to = "groups",values_to ="clust") %>%
      dplyr::mutate(groups = as.integer(groups)) %>%
      tidyr::nest(clusters = c(clust)) %>%
      dplyr::mutate(clusters = future_map(clusters
                                          , . %>%
                                            dplyr::mutate(!!ensym(siteCol) := siteNames
                                                          , cluster = numbers2words(clust)
                                                          , cluster = fct_reorder(cluster,clust)
                                                          )
                                          )
                    )
    
  }
  
  
  clustering_summarise <- function(clustDf,groupName="clust") {
    
    clust <- clustDf %>% dplyr::select(all_of(groupName))
    
    tab <- table(clust)
    
    tibble(nSites = nrow(clust)
           , nClusters = length(tab)
           , minClustSize = min(tab)
           , avClustSize = mean(tab)
           , maxClustSize = max(tab)
           )
    
  }
  
  
  clustering_explore <- function(clusters)
    
    
    
  

  # Make a cluster data frame (now prefer to use make_clusters)
  make_cluster_df <- function(rawClusters,siteIDsDf) {
    
    siteIDsDf %>%
      dplyr::mutate(id = row_number()) %>%
      dplyr::bind_cols(rawClusters %>%
                         dplyr::mutate(cluster = numbers2words(clust)
                                       , cluster = fct_reorder(cluster,clust)
                                       )
                       )
    
  }
  

  
  
  make_sil <- function(clustDf, distObj = datDist, clustID = "clust"){
    
    clustCol <- if(is.character(clustID)) which(names(clustDf)==clustID) else siteID
    
    silhouette(dplyr::pull(clustDf,clustCol),distObj)
    
  }
  
  # Turn an object of class silhouette into a data frame with one row per site
  make_sil_df <- function(clustDf,silObj) {
    
    clustDf %>%
      dplyr::bind_cols(tibble(neighbour = silObj[,2],sil_width = silObj[,3]))
    
  }
  
  # Calculate the within clusters sum of squares using only the floristic distances (datDist)
  calc_SS <- function(clustDf,clustDfMin,dist = datDist) {
    
    if(!exists("sqDist")) sqDist <- as.matrix(dist^2)
    
    clustDf %>%
      dplyr::mutate(id = row_number()) %>%
      dplyr::inner_join(clustDfMin) %>%
      dplyr::group_by(cluster) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      #dplyr::sample_n(2) %>% # TESTING
      dplyr::mutate(wss = map_dbl(data
                                  , ~sum(sqDist[.$id,.$id])/(2*nrow(.))
                                  )
                    ) %>%
      dplyr::select(-data)
    
  }
  
  
#-------Clustering - Indval--------
  
  cluster_indval <- function(clustDf,dfWithNames = datWide){
    
    # clustDf is usually siteID, clust, cluster
    # datWide is usually siteID * spp
    
    datForInd <- clustDf %>%
      dplyr::inner_join(dfWithNames) %>%
      .[,colSums(. != 0) > 0]
    
    labdsv::indval(datForInd[,names(datForInd) %in% names(datWide[,-1])]
                   ,datForInd$clust
                   )
    
  }
  
  
  # Indval result
  cluster_indval_df <- function(clustInd,clustDf){
    
    tibble(Taxa = names(clustInd$maxcls)
           , clust = clustInd$maxcls
           , indval = clustInd$indcls
           , pval = clustInd$pval
           ) %>%
      dplyr::inner_join(clustInd$relabu %>%
                          as_tibble(rownames = "Taxa") %>%
                          tidyr::gather(clust,abu,names(.)[names(.) %in% unique(clustDf$clust)]) %>%
                          dplyr::mutate(clust = as.numeric(clust)) %>%
                          dplyr::filter(abu > 0)
                       ) %>%
      dplyr::inner_join(clustInd$relfrq %>%
                         as_tibble(rownames = "Taxa") %>%
                         tidyr::gather(clust,frq,names(.)[names(.) %in% unique(clustDf$clust)]) %>%
                         dplyr::mutate(clust = as.numeric(clust)) %>%
                         dplyr::filter(frq > 0)
                       ) %>%
      dplyr::mutate(clust = as.numeric(clust)
                    , cluster = numbers2words(as.numeric(clust))
                    , cluster = fct_reorder(cluster,clust)
                    )
    
  }
  
  #------Clustering - Diagnostics----------
  
  
  diagnostic_plot <- function(diagnosticDf,labelDiagnostic = "diagnostic") {
    
    ggplot(diagnosticDf
           ,aes(groups
                , value
                , colour = combo
                , alpha = weight
                , label = groups
                , size = top
                )
           ) +
      geom_point() +
      geom_text_repel(data = diagnosticDf %>%
                        dplyr::filter(best)
                      , size = 2
                      , show.legend = FALSE
                      , box.padding = 1
                      , min.segment.length = 0
                      , colour = "black"
                      ) +
      facet_grid(as.formula(paste0(labelDiagnostic,"~method"))
                 , scales="free_y"
                 ,  labeller = label_wrap_gen(25,multi_line = TRUE)
                 ) +
      labs(colour = "Combination"
           , alpha = "Diagnostic used" #paste0("Top ",unique(diagnosticDf$topThresh)*100,"%")
           , title = paste0("Labels indicate top ",numbers2words(unique(diagnosticDf$bestThresh))," results")
           , size = paste0("Best ",(1-unique(diagnosticDf$topThresh))*100,"%")
           ) +
      scale_colour_viridis_c() +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      theme(strip.text.y = element_text(angle = 0))
    
  }
  
  
  diagnostic_df <- function(df
                            , useWeights = "weightClusters"
                            , diagnosticMinGroups = targetMinGroups
                            , summariseMethod = median
                            , topThresh = 2/3
                            , bestThresh = 5
                            , diagnosticDf = diagnostics
                            ) {
    
    df %>%
      dplyr::filter(groups > diagnosticMinGroups) %>%
      dplyr::select(method,groups,any_of(diagnostics$diagnostic)) %>%
      #select_if(~ !all(is.na(.))) %>%
      dplyr::group_by(method,groups) %>%
      dplyr::summarise(across(any_of(diagnostics$diagnostic),summariseMethod)) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_longer(any_of(diagnostics$diagnostic),names_to = "diagnostic", values_to = "value") %>%
      dplyr::left_join(diagnostics) %>%
      dplyr::group_by(diagnostic,diagDefinition) %>%
      dplyr::mutate(scale = if_else(highGood
                                    ,scales::rescale(value,to=c(0,1))
                                    ,scales::rescale(desc(value),to=c(0,1))
                                    )
                    ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(comboInit = scale*!!ensym(useWeights)) %>%
      dplyr::group_by(method,groups) %>%
      dplyr::mutate(combo = mean(comboInit)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(combo = scales::rescale(combo,to=c(0,1))) %>%
      dplyr::mutate(topThresh = topThresh
                    , bestThresh = bestThresh
                    , top = combo >= topThresh
                    , top = if_else(is.na(top),FALSE,top)
                    , best = combo >= sort(unique(.$combo),TRUE)[bestThresh]
                    , best = if_else(is.na(best),FALSE,best)
                    , diagnostic = factor(diagnostic, levels = levels(diagnostics$diagnostic))
                    , weight = !!ensym(useWeights)
                    )
    
  }
  
#-------


  dist_to_df <- function(inDist,patchNames) {
    
    # https://stackoverflow.com/questions/23474729/convert-object-of-class-dist-into-data-frame-in-r/23475065
    
    if (class(inDist) != "dist") stop("wrong input type")
    A <- attr(inDist, "Size")
    B <- patchNames
    if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
    if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
    data.frame(
      row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
      col = rep(B[-length(B)], (length(B)-1):1),
      value = as.vector(inDist))
    
  }


  


# Convert deg°min'sec to decimal degrees
  deg_min_sec_to_dec_deg <- function(df,col,sep = "\u00B0|'", into = c("lat","long")) {
    
    dfNames <- names(df)
    
    df %>%
      tidyr::separate(!!ensym(col), into = c("deg","min","sec"), sep = paste0(sep)) %>%
      dplyr::mutate(sec = parse_number(sec)) %>%
      dplyr::mutate(across(all_of(c("deg","min","sec")),as.numeric)) %>%
      dplyr::mutate(!!ensym(into) := deg + min/60 + sec/60^2) %>%
      dplyr::select(all_of(dfNames),!!ensym(into))
    
  }
  
# get a bib entry from a DOI
  get_bib <- function(DOI,outFile = NULL){
    
    # https://stackoverflow.com/questions/57340204/r-convert-list-of-dois-to-bibtex
    
    h <- curl::new_handle()
    curl::handle_setheaders(h, "accept" = "application/x-bibtex")
    
    get_the_bib <- function(doi) {
      
      try(
        curl::curl(doi,handle=h) %>%
        readLines(warn = FALSE) %>%
        {if(is.character(outFile)) write(.,file=outFile,append = TRUE) else (.)}
      )
      
    }
    
    DOI %>%
      gsub("https://doi.org/","",.,fixed = TRUE) %>%
      paste0("https://doi.org/",.) %>%
      map(get_the_bib)
    
  }


  # Git add.
  gitadd <- function(dir = getwd()){
    cmd_list <- list(
      cmd1 = tolower(substr(dir,1,2)),
      cmd2 = paste("cd",dir),
      cmd3 = "git add --all"
    )
    cmd <- paste(unlist(cmd_list),collapse = " & ")
    shell(cmd)
  }
  
  # Git commit.
  gitcommit <- function(msg = "commit from Rstudio", dir = getwd()){
    cmd = sprintf("git commit -m\"%s\"",msg)
    shell(cmd)
  }
  
  # Git push.
  gitpush <- function(dir = getwd()){
    cmd_list <- list(
      cmd1 = tolower(substr(dir,1,2)),
      cmd2 = paste("cd",dir),
      cmd3 = "git push"
    )
    cmd <- paste(unlist(cmd_list),collapse = " & ")
    shell(cmd)
  }


# Site a package in rmarkdown
  # assumes these have been run and that 'packages' contains all packages to cite
    # knnitr::write_bib(packages,"packageCitations.bib")
    # refs <- bib2df::bib2df("packageCitations.bib")

  cite_package <- function(package,brack = TRUE,startText = "", endText = "") {
    
    thisRef <- refs %>%
      dplyr::filter(grepl(paste0("-",package),BIBTEXKEY) | grepl(paste0("^",package),BIBTEXKEY)) %>%
      dplyr::pull(BIBTEXKEY)
    
    starts <- if(brack) paste0("[",startText,"@") else paste0(startText,"@")
    ends <- if(brack) paste0(endText,"]") else endText
    
    if(length(thisRef) > 1) {
      
      paste0(starts,paste0(thisRef,collapse = "; @"),ends)
      
    } else {
      
      paste0(starts,"R-",package,ends)
      
    }
    
  }
  
# make package bibliography, including tweaks for known package issues.
  fix_bib <- function(bibFile, makeKey = FALSE, isPackageBib = FALSE) {
    
    inRefs <- bib2df::bib2df(bibFile)
    
    namesInRefs <- colnames(inRefs) %>%
      grep("\\.\\d+$",.,value = TRUE,invert = TRUE) %>%
      `[`(1:28) %>%
      c(.,"COPYRIGHT")
    
    refs <- inRefs %>%
      {if(isPackageBib) (.) %>% dplyr::mutate(package = gsub("R-|\\d{4}","",BIBTEXKEY)) else (.)} %>%
      tidytext::unnest_tokens("titleWords"
                              ,TITLE
                              ,token = "regex"
                              ,pattern = " "
                              ,to_lower = FALSE
                              #,strip_punct = FALSE
                              ,collapse = FALSE
                              ) %>%
      dplyr::mutate(titleWords = gsub("\\{|\\}","",titleWords)
                    , isCap = grepl(paste0(LETTERS,collapse="|"),titleWords)
                    , titleWords = if_else(isCap,paste0("{",titleWords,"}"),titleWords)
                    ) %>%
      tidyr::nest(data = c(titleWords,isCap)) %>%
      dplyr::mutate(TITLE = map_chr(data,. %>% dplyr::pull(titleWords) %>% paste0(collapse = " "))
                    , AUTHOR = map(AUTHOR,~gsub("Microsoft Corporation","{Microsoft Corporation}",.))
                    , AUTHOR = map(AUTHOR,~gsub("Fortran original by |R port by ","",.))
                    , AUTHOR = map(AUTHOR, ~gsub("with contributions by","and",.))
                    , AUTHOR = map(AUTHOR, ~gsub("Â "," ",.))
                    , YEAR = substr(YEAR,1,4)
                    ) %>%
      {if(makeKey) (.) %>%
          dplyr::mutate(BIBTEXKEY = map2_chr(AUTHOR
                                       ,YEAR
                                       ,~paste0(toupper(gsub("[[:punct:]]|\\s","",.x[[1]]))
                                                , .y
                                                )
                                       )
                        ) else (.)} %>%
      {if(isPackageBib) (.) %>% dplyr::mutate(TITLE = map2_chr(package,TITLE,~gsub(.x,paste0("{",.x,"}"),.y))) else (.)} %>%
      dplyr::select(any_of(namesInRefs)) %>%
      dplyr::filter(!grepl("MEDIA SCREEN AND",CATEGORY))
    
    bib2df::df2bib(refs,bibFile)
    
    return(refs)
    
  }


# Are the values within a column unique
  col_is_unique <- function(df,col = "SiteID") {
    
    notUnique <- df %>%
      dplyr::select(grep("^n$",names(.),value = TRUE,invert = TRUE)) %>%
      dplyr::count(!!ensym(col)) %>%
      dplyr::filter(n > 1)
    
    print(paste0("there are ",nrow(notUnique)," ",col,"(s) that are not unique: ",vec_to_sentence(notUnique[,1])))
    
  }


# Unscale scaled data
unscale_data <- function(scaledData) {
  
  scaledData*attr(scaledData,"scaled:scale")+attr(scaledData,"scaled:center")
  
}


# From https://github.com/mtennekes/tmap/issues/255 - who got it from:
#http://stackoverflow.com/questions/20241065/r-barplot-wrapping-long-text-labels

# Core wrapping function
  wrap.it <- function(x, len) { 
    
    sapply(x, function(y) paste(strwrap(y, len), 
                                collapse = "\n"), 
           USE.NAMES = FALSE)
  }

# Call this function with a list or vector
  wrap.labels <- function(x, len) {
    
    if (is.list(x))
    {
      lapply(x, wrap.it, len)
    } else {
      wrap.it(x, len)
    }
  }



# From https://gist.github.com/danlwarren/271288d5bab45d2da549

# Function to rarefy point data in any number of dimensions.  The goal here is to 
# take a large data set and reduce it in size in such a way as to approximately maximize the 
# difference between points.  For instance, if you have 2000 points but suspect a lot of 
# spatial autocorrelation between them, you can pass in your data frame, the names (or indices)
# of the lat/lon columns, and the number 200, and you get back 200 points from your original data 
# set that are chosen to be as different from each other as possible given a randomly chosen
# starting point

# Input is:
#
# x, a data frame containing the columns to be used to calculate distances along with whatever other data you need
# cols, a vector of column names or indices to use for calculating distances
# npoints, the number of rarefied points to spit out
#
# e.g., thin.max(my.data, c("latitude", "longitude"), 200)
  
  
  thin_max <- function(x, cols, npoints){
    #Create empty vector for output
    inds <- vector(mode="numeric")
    
    #Create distance matrix
    this.dist <- as.matrix(dist(x[,cols], upper=TRUE))
    
    #Draw first index at random
    inds <- c(inds, as.integer(runif(1, 1, length(this.dist[,1]))))
    
    #Get second index from maximally distant point from first one
    #Necessary because apply needs at least two columns or it'll barf
    #in the next bit
    inds <- c(inds, which.max(this.dist[,inds]))
    
    while(length(inds) < npoints){
      #For each point, find its distance to the closest point that's already been selected
      min.dists <- apply(this.dist[,inds], 1, min)
      
      #Select the point that is furthest from everything we've already selected
      this.ind <- which.max(min.dists)
      
      #Get rid of ties, if they exist
      if(length(this.ind) > 1){
        print("Breaking tie...")
        this.ind <- this.ind[1]
      }
      inds <- c(inds, this.ind)
    }
    
    return(x[inds,])
  }

# A function to run random forest over a df with first column 'cluster' and other columns explanatory

  rf_mod_fold <- function(envClust
                          , clustCol = "cluster"
                          , envCols = names(patchesEnvSelect)[-1]
                          , idCol = "cell"
                          , doFolds = folds
                          , outFile
                          , saveModel = FALSE
                          , saveImp = FALSE
                          , ...
                          ){
    
    idCol <- if(is.numeric(idCol)) names(envClust)[idCol] else idCol
    
    clustCol <- if(is.numeric(clustCol)) names(envClust)[clustCol] else clustCol
    
    envCols <- if(is.numeric(envCols)) names(envClust)[envCols] else envCols
    
    envClust <- envClust %>%
      dplyr::mutate(fold = sample(1:doFolds,nrow(.),replace=TRUE,rep(1/doFolds,doFolds)))
    
    folds <- 1:doFolds
    
    fold_rf_mod <- function(fold) {
      
      outFile <- gsub("_conf",paste0("_fold",fold,"_conf"),outFile)
      
      if(!file.exists(outFile)) {
        
        if(doFolds > 1) {
          
          train <- envClust[envClust$fold != fold,which(names(envClust) %in% c(clustCol,envCols))] %>%
            dplyr::mutate(cluster = factor(cluster))
          
          test <- envClust[envClust$fold == fold,which(names(envClust) %in% c(idCol,clustCol,envCols))] %>%
            dplyr::mutate(cluster = factor(cluster, levels = levels(train$cluster)))
          
          rfMod <- randomForest(x = train[,envCols]
                                , y = train[,clustCol] %>% dplyr::pull(!!ensym(clustCol))
                                , ntree = 500
                                , importance = saveImp
                                )
        
          rfPred <- test %>%
            dplyr::select(!!ensym(idCol),!!ensym(clustCol)) %>%
            dplyr::bind_cols(predict(rfMod
                                     , newdata = test[,envCols]
                                     ) %>%
                               tibble::enframe(name = NULL, value = "predCluster")
                             ) %>%
            dplyr::bind_cols(predict(rfMod
                                     , newdata = test[,envCols]
                                     , type = "prob"
                                     ) %>%
                               as_tibble()
                             )
          
        } else {
          
          rfMod <- randomForest::randomForest(x = envClust[,which(names(envClust) %in% c(envCols))]
                                              , y = envClust %>% dplyr::pull(!!ensym(clustCol))
                                              , ntree = 500
                                              , importance = saveImp
                                              )
          
          rfPred <- envClust[,c(idCol,clustCol)] %>%
            dplyr::bind_cols(predict(rfMod) %>%
                               tibble::enframe(name = NULL, value = "predCluster")
                             ) %>%
            dplyr::bind_cols(predict(rfMod
                                     , type = "prob"
                                     ) %>%
                               as_tibble()
                             )
          
        }
        
        feather::write_feather(rfPred,outFile)
        
        if(saveImp) {feather::write_feather(as_tibble(rfMod$importance, rownames = "att"),gsub("_rfPred","_rfImp",outFile))}
        
        if(saveModel) {feather::write_feather(rfMod,gsub("_rfPred","",outFile))}
        
      }
      
    }
    
    map(folds,fold_rf_mod)
    
    }
  
  rf_mod <- function(envClust
                          , clustCol
                          , envCols
                          , idCol
                          , outFile
                          , saveModel = FALSE
                          , saveImp = FALSE
                          , ...
                     ){
    
    idCol <- if(is.numeric(idCol)) names(envClust)[idCol] else idCol
    
    clustCol <- if(is.numeric(clustCol)) names(envClust)[clustCol] else clustCol
    
    envCols <- if(is.numeric(envCols)) names(envClust)[envCols] else envCols
    
    rfMod <- randomForest::randomForest(x = envClust[,which(names(envClust) %in% c(envCols))]
                                        , y = envClust %>% dplyr::pull(!!ensym(clustCol))
                                        , ntree = 500
                                        , importance = saveImp
                                        )
    
    rfPred <- envClust[,c(idCol,clustCol)] %>%
      dplyr::bind_cols(predict(rfMod) %>%
                         tibble::enframe(name = NULL, value = "predCluster")
                       ) %>%
      dplyr::bind_cols(predict(rfMod
                               , type = "prob"
                               ) %>%
                         as_tibble()
                       )
    
    feather::write_feather(rfPred,outFile)
    
    if(saveImp) {feather::write_feather(as_tibble(rfMod$importance, rownames = "att"),gsub("_rfPred","_rfImp",outFile))}
    
    if(saveModel) {feather::write_feather(rfMod,gsub("_rfPred","",outFile))}
    
  }


# A function to run random forest using tidymodels dialogue
  
  run_rf <- function(datTrain,modRecipe) {
    
    randomForest::randomForest(as.formula(modRecipe)
                               , data = datTrain
                               , ntree = 500
                               )
    
  }
    
    
# Function to get data out of 32 bit MS Access from 64 bit R
# see https://stackoverflow.com/questions/13070706/how-to-connect-r-with-access-database-in-64-bit-window

access_query_32 <- function(db_path, db_table = "qryData_RM", table_out = "data_access") {
  
  library(svSocket)
  
  # variables to make values uniform
  sock_port <- 8642L
  sock_con <- "sv_con"
  ODBC_con <- "a32_con"
  
  if (file.exists(db_path)) {
    
    # build ODBC string
    ODBC_str <- local({
      s <- list()
      s$path <- paste0("DBQ=", gsub("(/|\\\\)+", "/", path.expand(db_path)))
      s$driver <- "Driver={Microsoft Access Driver (*.mdb, *.accdb)}"
      s$threads <- "Threads=4"
      s$buffer <- "MaxBufferSize=4096"
      s$timeout <- "PageTimeout=5"
      paste(s, collapse=";")
    })
    
    # start socket server to transfer data to 32 bit session
    startSocketServer(port=sock_port, server.name="access_query_32", local=TRUE)
    
    # build expression to pass to 32 bit R session
    expr <- "library(svSocket)"
    expr <- c(expr, "library(RODBC)")
    expr <- c(expr, sprintf("%s <- odbcDriverConnect('%s')", ODBC_con, ODBC_str))
    expr <- c(expr, sprintf("if('%1$s' %%in%% sqlTables(%2$s)$TABLE_NAME) {%1$s <- sqlFetch(%2$s, '%1$s')} else {%1$s <- 'table %1$s not found'}", db_table, ODBC_con))
    expr <- c(expr, sprintf("%s <- socketConnection(port=%i)", sock_con, sock_port))
    expr <- c(expr, sprintf("evalServer(%s, %s, %s)", sock_con, table_out, db_table))
    expr <- c(expr, "odbcCloseAll()")
    expr <- c(expr, sprintf("close(%s)", sock_con))
    expr <- paste(expr, collapse=";")
    
    # launch 32 bit R session and run expressions
    prog <- file.path(R.home(), "bin", "i386", "Rscript.exe")
    system2(prog, args=c("-e", shQuote(expr)), stdout=NULL, wait=TRUE, invisible=TRUE)
    
    # stop socket server
    stopSocketServer(port=sock_port)
    
    # display table fields
    message("retrieved: ", table_out, " - ", paste(colnames(get(table_out)), collapse=", "))
  } else {
    warning("database not found: ", db_path)
  }
}

# Function to create grid of continuous values - usually for prediction

  create_grid_cont <- function(df,vars,seqLength){
    
    df %>%
      dplyr::summarise_at(vars, min) %>%
      dplyr::mutate(type = "min") %>%
      dplyr::bind_rows(df %>%
                         dplyr::summarise_at(vars, max) %>%
                         dplyr::mutate(type = "max")
                       ) %>%
      tidyr::gather(variable,value,1:length(vars)) %>%
      tidyr::spread(type,value) %>%
      dplyr::mutate(values = map2(min,max,~seq(.x,.y,(.y-.x)/(seqLength-1)))) %>%
      dplyr::select(variable,values) %>%
      tidyr::unnest() %>%
      split(.$variable) %>%
      lapply(function(x) x %>% dplyr::select(2) %>% unlist()) %>%
      expand.grid()
    
  }

# Function written by Andrew Bevan, found on R-sig-Geo, and modified by Pascal Title
ah2sp <- function(x, increment=360, rnd=10, proj4string=CRS(as.character(NA)),tol=1e-4) {
  if (class(x) != "ahull"){
    stop("x needs to be an ahull class object")
  }
  # Extract the edges from the ahull object as a dataframe
  xdf <- as.data.frame(x$arcs)
  
  #correct for possible arc order strangeness (Pascal Title addition 29 Nov 2013)
  k <- 1
  xdf <- cbind(xdf, flip = rep(FALSE, nrow(xdf)))
  repeat{
    if (is.na(xdf[k+1, 'end1'])) {
      break
    }
    #cat(k, '\n')
    if (xdf[k,'end2'] == xdf[k+1,'end1']) {
      #cat('case 1\n')
      k <- k + 1
    } else if (xdf[k,'end2'] != xdf[k+1,'end1'] & !xdf[k,'end2'] %in% xdf$end1[k+1:nrow(xdf)] & !xdf[k,'end2'] %in% xdf$end2[k+1:nrow(xdf)]) {
      #cat('case 2\n')
      k <- k + 1
    } else if (xdf[k,'end2'] != xdf[k+1,'end1'] & xdf[k,'end2'] %in% xdf$end1[k+1:nrow(xdf)] & !xdf[k,'end2'] %in% xdf$end2[k+1:nrow(xdf)]) {
      #cat('case 3\n')
      m <- which(xdf$end1[k+1:nrow(xdf)] == xdf[k,'end2']) + k
      xdf <- rbind(xdf[1:k,],xdf[m,],xdf[setdiff((k+1):nrow(xdf),m),])
    } else if (xdf[k,'end2'] != xdf[k+1,'end1'] & !xdf[k,'end2'] %in% xdf$end1[k+1:nrow(xdf)] & xdf[k,'end2'] %in% xdf$end2[k+1:nrow(xdf)]) {
      #cat('case 4\n')
      m <- which(xdf$end2[k+1:nrow(xdf)] == xdf[k,'end2']) + k
      tmp1 <- xdf[m,'end1']
      tmp2 <- xdf[m,'end2']
      xdf[m,'end1'] <- tmp2
      xdf[m,'end2'] <- tmp1
      xdf[m,'flip'] <- TRUE
      xdf <- rbind(xdf[1:k,], xdf[m,], xdf[setdiff((k+1):nrow(xdf), m),])
    } else {
      k <- k + 1
    }
  }	
  
  
  # Remove all cases where the coordinates are all the same      
  xdf <- subset(xdf, xdf$r > 0)
  res <- NULL
  if (nrow(xdf) > 0) {
    # Convert each arc to a line segment
    linesj <- list()
    prevx <- NULL
    prevy <- NULL
    j <- 1
    for(i in 1:nrow(xdf)) {
      rowi <- xdf[i,]
      v <- c(rowi$v.x, rowi$v.y)
      theta <- rowi$theta
      r <- rowi$r
      cc <- c(rowi$c1, rowi$c2)
      # Arcs need to be redefined as strings of points. Work out the number of points to allocate in this arc segment.
      ipoints <- 2 + round(increment * (rowi$theta / 2), 0)
      # Calculate coordinates from arc() description for ipoints along the arc.
      angles <- alphahull::anglesArc(v, theta)
      if (rowi['flip'] == TRUE){ angles <- rev(angles) }
      seqang <- seq(angles[1], angles[2], length = ipoints)
      x <- round(cc[1] + r * cos(seqang),rnd)
      y <- round(cc[2] + r * sin(seqang),rnd)
      # Check for line segments that should be joined up and combine their coordinates
      if (is.null(prevx)) {
        prevx <- x
        prevy <- y
        # added numerical precision fix (Pascal Title Dec 9 2013)
      } else if ((x[1] == round(prevx[length(prevx)],rnd) | abs(x[1] - prevx[length(prevx)]) < tol) && (y[1] == round(prevy[length(prevy)],rnd) | abs(y[1] - prevy[length(prevy)]) < tol)) {
        if (i == nrow(xdf)){
          #We have got to the end of the dataset
          prevx <- append(prevx ,x[2:ipoints])
          prevy <- append(prevy, y[2:ipoints])
          prevx[length(prevx)] <- prevx[1]
          prevy[length(prevy)] <- prevy[1]
          coordsj <- cbind(prevx,prevy)
          colnames(coordsj) <- NULL
          # Build as Line and then Lines class
          linej <- Line(coordsj)
          linesj[[j]] <- Lines(linej, ID = as.character(j))
        } else {
          prevx <- append(prevx, x[2:ipoints])
          prevy <- append(prevy, y[2:ipoints])
        }
      } else {
        # We have got to the end of a set of lines, and there are several such sets, so convert the whole of this one to a line segment and reset.
        prevx[length(prevx)] <- prevx[1]
        prevy[length(prevy)] <- prevy[1]
        coordsj <- cbind(prevx,prevy)
        colnames(coordsj)<-NULL
        # Build as Line and then Lines class
        linej <- Line(coordsj)
        linesj[[j]] <- Lines(linej, ID = as.character(j))
        j <- j + 1
        prevx <- NULL
        prevy <- NULL
      }
    }
    
    #Drop lines that will not produce adequate polygons (Pascal Title addition 9 Dec 2013)
    badLines <- vector()
    for (i in 1:length(linesj)){
      if (nrow(linesj[[i]]@Lines[[1]]@coords) < 4){
        badLines <- c(badLines,i)
      }
    }
    if (length(badLines) > 0){linesj <- linesj[-badLines]}
    
    # Promote to SpatialLines
    lspl <- SpatialLines(linesj)
    # Convert lines to polygons
    # Pull out Lines slot and check which lines have start and end points that are the same
    lns <- slot(lspl, "lines")
    polys <- sapply(lns, function(x) { 
      crds <- slot(slot(x, "Lines")[[1]], "coords")
      identical(crds[1, ], crds[nrow(crds), ])
    }) 
    # Select those that do and convert to SpatialPolygons
    polyssl <- lspl[polys]
    list_of_Lines <- slot(polyssl, "lines")
    sppolys <- SpatialPolygons(list(Polygons(lapply(list_of_Lines, function(x) { Polygon(slot(slot(x, "Lines")[[1]], "coords")) }), ID = "1")), proj4string=proj4string)
    # Create a set of ids in a dataframe, then promote to SpatialPolygonsDataFrame
    hid <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "ID"))
    areas <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "area"))
    df <- data.frame(hid,areas)
    names(df) <- c("HID","Area")
    rownames(df) <- df$HID
    res <- SpatialPolygonsDataFrame(sppolys, data=df)
    res <- res[which(res@data$Area > 0),]
  }  
  return(res)
}


# glm 'chi-square'

  contingency_glm <- function(cont){
    
    var1 <- names(cont)[1]
    var2 <- names(cont)[2]
    
    contingency <- cont %>%
      dplyr::mutate(var1 = factor(!!ensym(var1))
                    , var2 = factor(!!ensym(var2))
                    ) %>%
      dplyr::group_by(var2) %>%
      dplyr::filter(success > 3) %>%
      dplyr::mutate(levels = n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(levels == max(levels)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(!!var1 := if_else(is.na(!!ensym(var1)),var1,!!ensym(var1))
                    , !!var2 := if_else(is.na(!!ensym(var2)),var2,!!ensym(var2))
                    , var1 = fct_inorder(factor(as.character(var1)))
                    , var2 = fct_reorder(as.character(var2), success)
                    , var1No = as.factor(as.numeric(var1))
                    , var2No = as.factor(as.numeric(var2))
                    , trials = success + failure
                    ) %>%
      replace(is.na(.), 0)
    
    
    mod <- stan_glm(cbind(success,failure) ~ var1*var2
                    , data = contingency
                    , family = binomial()
                    #, iter = 5000
                    )
    
    modFit <- pp_check(mod)
    
    modRhat <- plot(mod, "rhat_hist")
    
    modTrace <- stan_trace(mod)
    
    mod2d <- pp_check(mod, plotfun = "stat_2d", stat = c("mean", "sd"))
    
    # Use the model to predict results over variables of interest
    modPred <- contingency %>%
      dplyr::group_by(var1,var2) %>%
      dplyr::summarise(success = 0, failure = 100) %>% # use 100 trials to give results as percentages
      dplyr::ungroup() %>%
      dplyr::mutate(col = row.names(.)) %>%
      dplyr::left_join(as_tibble(posterior_predict(mod
                                                   , newdata = .
                                                   )
                                 ) %>%
                         tibble::rownames_to_column(var = "row") %>%
                         tidyr::gather(col,value,2:ncol(.))
                       ) %>%
      dplyr::left_join(contingency %>%
                         dplyr::select(!!ensym(var1),!!ensym(var2),var1,var2) %>%
                         unique()
                       )
    
    # summarise the results
    modRes <- as_tibble(modPred) %>%
      dplyr::group_by(!!ensym(var1),!!ensym(var2)) %>%
      dplyr::summarise(n = n()
                       , nCheck = nrow(as_tibble(mod))
                       , modMedian = quantile(value,0.5)
                       , modMean = mean(value)
                       , modci90lo = quantile(value, 0.025)
                       , modci90up = quantile(value, 0.975)
                       , ci = modci90up-modci90lo
                       , text = paste0(round(modMedian,2)," (",round(modci90lo,2)," to ",round(modci90up,2),")")
                       ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate_if(is.numeric,round,2)
    
    
    modPlotRidges <- ggplot(modPred, aes(value,!!ensym(var1),fill=!!ensym(var1))) +
      ggridges::geom_density_ridges(alpha = 0.5) +
      facet_wrap(~get("var2"),scales="free") +
      scale_fill_viridis_d()
    
    modDiffVar1 <- modPred %>%
      (function(x) x %>% dplyr::left_join(x %>% dplyr::select(var1,var2,row,value) %>% dplyr::rename(var1b = var1, value_2 = value))) %>%
      dplyr::filter(var1 != var1b) %>%
      dplyr::mutate(comparison = map2_chr(var1,var1b,~paste(sort(c(.x,.y))[1],sort(c(.x,.y))[2]))) %>%
      dplyr::group_by(comparison,var2,row) %>%
      dplyr::slice(1) %>% # from memory this is part of a trick to remove duplicates... comparison field is key to this?
      dplyr::ungroup() %>%
      dplyr::mutate(diff = value-value_2
                    , maxEst = map2_dbl(abs(value),abs(value_2),max)
                    , test01 = map2_lgl(value,value_2,~all(.x==0,.y==0))
                    , perDiff = if_else(!test01,100*diff/maxEst,0)
                    )
    
    modDiffRes <- modDiffVar1 %>%
      dplyr::group_by(var1,var2,var1b) %>%
      dplyr::summarise(n = n()
                       , nCheck = nrow(as_tibble(mod))
                       , meanDiff = median(diff)
                       , value = mean(perDiff)
                       ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(alpha = 1.5*abs(value-50)
                    , colour = if_else(value > 0,"More likely","Less likely")
                    , text = paste0(round(abs(value),0),"%")
                    , longText = paste0(var2, " was ", text, " ", tolower(colour), " to occur in ", var1, " reviews than ",var1b," reviews")
                    , var1Comparison = paste0(var1," vs ",var1b)
                    , var1Comparison = fct_relevel(var1Comparison, grep("SA",unique(var1Comparison),value=TRUE))
                    ) %>%
      dplyr::mutate_if(is.numeric,round,2) %>%
      dplyr::arrange(var2,var1,var1b)
    
    modPlotRidgesDiff <- ggplot(modDiffVar1, aes(perDiff,paste0(get("var1")," vs ",get("var1b")))) +
      ggridges::geom_density_ridges(alpha = 0.5) +
      geom_vline(aes(xintercept = 0)) +
      facet_wrap(~get("var2"),scales="free") +
      scale_fill_viridis_d() +
      labs(y = "Comparison"
           , x = "Difference"
           )
    
    modPlotDiff <- ggplot(modDiffRes,aes(var1Comparison,var2,fill=colour,alpha=abs(value),label=text)) +
      geom_tile() +
      geom_text(size=3) +
      scale_fill_viridis_d() +
      scale_alpha_continuous(guide = FALSE
                             , range = c(0,0.5)
                             ) +
      scale_x_discrete(limits = rev(levels(var2))) +
      labs(fill = "Likelihood"
           , x = "Comparison"
           , y = get("var2")
           )
    
     mget(ls(pattern = "mod"))
    
  }

  # chi-square
  chi_square <- function(cont) {
    
    var1 <- names(cont)[1]
    var2 <- names(cont)[2]
    
    contingency <- cont %>%
      dplyr::mutate(var1 = factor(!!ensym(var1))
                    , var2 = factor(!!ensym(var2))
                    ) %>%
      tidyr::complete(var1,var2) %>%
      dplyr::mutate(!!var1 := if_else(is.na(!!ensym(var1)),var1,!!ensym(var1))
                    , !!var2 := if_else(is.na(!!ensym(var2)),var2,!!ensym(var2))
                    , var1 = fct_inorder(factor(as.character(var1)))
                    , var2 = fct_inorder(factor(as.character(var2)))
                    , var1No = as.factor(as.numeric(var1))
                    , var2No = as.factor(as.numeric(var2))
                    ) %>%
      replace(is.na(.), 0)
    
    chSq <- contingency %>%
      dplyr::select(var1No,var2No,n) %>%
      tidyr::spread(var2No,n,drop=TRUE) %>%
      as.data.frame %>%
      tibble::column_to_rownames(names(.)[1]) %>%
      chisq.test()
    
    chSqResidual <- chSq$residuals %>%
      data.frame %>%
      tibble::rownames_to_column("var1No") %>%
      tidyr::gather("var2No","residual",2:ncol(.)) %>%
      dplyr::mutate(var2No = gsub("X","",var2No))
    
    chSqVis <- data.frame(100*chSq$residuals^2/chSq$statistic) %>%
      data.frame %>%
      tibble::rownames_to_column("var1No") %>%
      tidyr::gather("var2No","contribution",2:ncol(.))%>%
      dplyr::mutate(var2No = gsub("X","",var2No)) %>%
      as_tibble() %>%
      dplyr::left_join(chSqResidual) %>%
      dplyr::left_join(contingency) %>%
      dplyr::mutate(per = 100*contribution/sum(contribution)
                    , text = paste0("n:",n,"\n",round(per,1),"%")
                    , direction = if_else(residual<0
                                          ,"less than expected"
                                          , if_else(residual>0
                                                    ,"more than expected"
                                                    ,"as expected"
                                          )
                    )
                    , label = paste0(var2, " occurs ", direction, " in ", var1)
      ) %>%
      dplyr::select(!!var1,!!var2,contribution,residual,n,per,text,label,direction)
    
    chSqPlot <- ggplot(chSqVis, aes(!!ensym(var1), fct_rev(!!ensym(var2)), fill = direction, alpha = contribution, label = text)) +
      geom_tile() +
      geom_text(size = 2) +
      guides(alpha = FALSE) +
      theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
      labs(subtitle = "Percentages are the percent contribution to overall chi-squared value"
           , y = var2
           , x = var1
      ) +
      scale_fill_viridis_d()
    
    chSqText <- paste0("(Chi-squared = ",round(chSq$statistic,1), ", df = ",chSq$parameter,", p <= ",round(chSq$p.value,4),")")
    
    doChSqPlot <- chSq$p.value<0.05
    
    chSqRes <- list(chSq=chSq,chSqVis=chSqVis,chSqPlot=chSqPlot,chSqText=chSqText,doChSqPlot=doChSqPlot)
    
  }
  
# function to read in previously saved rds
  
  read_rds_file <- function(fileName) if(file.exists(fileName)) read_rds(fileName) else NULL
  
  
# Find peaks in a vector
  
  which_peaks <- function(x,partial=TRUE,decreasing=FALSE){
    if (decreasing){
      if (partial){
        which(diff(c(FALSE,diff(x)>0,TRUE))>0)
      }else {
        which(diff(diff(x)>0)>0)+1
      }
    } else {
      if (partial){
        which(diff(c(TRUE,diff(x)>=0,FALSE))<0)
      }else {
        which(diff(diff(x)>=0)<0)+1
      }
    }
  }
  
  
  
# Set column widths for kable tables
  
  html_table_width <- function(kable_output, width){
    width_html <- paste0(paste0('<col width="', width, '">'), collapse = "\n")
    sub("<table>", paste0("<table>\n", width_html), kable_output)
  }
  
  
  
  
# Create a colour palette for n groups

  col_pal <-  function(n) {
    if (n <= 8) {
      RColorBrewer::brewer.pal(n, "Set2")
    } else {
      hcl(h=seq(0,(n-1)/(n),length=n)*360,c=100,l=65,fixup=TRUE)
    }
  }
  

# turn a vector into a comma separated list of values with a penultimate 'and'
  vec_to_sentence <- function(x,sep=",",end="and") {
    
    x[!is.na(x)] %>%
      paste(collapse = "JOINSRUS") %>%
      (function(x) if(sep == ";") {
        
        stringi::stri_replace_last_regex(x,"JOINSRUS", paste0(sep," and ")) %>%
          str_replace_all("JOINSRUS",paste0(sep," "))
        
      } else {
        
        stringi::stri_replace_last_regex(x,"JOINSRUS",paste0(" ",end," ")) %>%
          str_replace_all("JOINSRUS",paste0(sep," "))
        
      }
      )
    
  }

# https://github.com/ateucher/useful_code/blob/master/R/numbers2words.r

  numbers2words <- function(x){
    ## Function by John Fox found here:
    ## http://tolstoy.newcastle.edu.au/R/help/05/04/2715.html
    ## Tweaks by AJH to add commas and "and"
    helper <- function(x){

      digits <- rev(strsplit(as.character(x), "")[[1]])
      nDigits <- length(digits)
      if (nDigits == 1) as.vector(ones[digits])
      else if (nDigits == 2)
        if (x <= 19) as.vector(teens[digits[1]])
      else trim(paste(tens[digits[2]],
                      Recall(as.numeric(digits[1]))))
      else if (nDigits == 3) trim(paste(ones[digits[3]], "hundred and",
                                        Recall(makeNumber(digits[2:1]))))
      else {
        nSuffix <- ((nDigits + 2) %/% 3) - 1
        if (nSuffix > length(suffixes)) stop(paste(x, "is too large!"))
        trim(paste(Recall(makeNumber(digits[
          nDigits:(3*nSuffix + 1)])),
          suffixes[nSuffix],"," ,
          Recall(makeNumber(digits[(3*nSuffix):1]))))
      }
    }
    trim <- function(text){
      #Tidy leading/trailing whitespace, space before comma
      text=gsub("^\ ", "", gsub("\ *$", "", gsub("\ ,",",",text)))
      #Clear any trailing " and"
      text=gsub(" and$","",text)
      #Clear any trailing comma
      gsub("\ *,$","",text)
    }
    makeNumber <- function(...) as.numeric(paste(..., collapse=""))
    #Disable scientific notation
    opts <- options(scipen=100)
    on.exit(options(opts))
    ones <- c("", "one", "two", "three", "four", "five", "six", "seven",
              "eight", "nine")
    names(ones) <- 0:9
    teens <- c("ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen",
               "sixteen", "seventeen", "eighteen", "nineteen")
    names(teens) <- 0:9
    tens <- c("twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty",
              "ninety")
    names(tens) <- 2:9
    x <- round(x)
    suffixes <- c("thousand", "million", "billion", "trillion")
    if (length(x) > 1) return(trim(sapply(x, helper)))
    res <- helper(x)
    #res <- gsub(" ","",res)
    return(res)
  }
  

# Fix html widget when not displayed
  widgetFix <- function(inputFile,outputFile){
    a = readLines(inputFile)
    output = paste(a,collapse="\n")
    output = gsub(">\n\n</div>","></div>",output)
    writeLines(output,outputFile)
    invisible(NULL)
  }
  

# Generate jenks breaks
# http://cainarchaeology.weebly.com/r-function-for-plotting-jenks-natural-breaks-classification.html
  
  plotJenks <- function(data, n=3, brks.cex=0.70, top.margin=10, dist=5){
    df <- data.frame(sorted.values=sort(data, decreasing=TRUE))
    Jclassif <- classIntervals(df$sorted.values, n, style = "jenks") #requires the 'classInt' package
    test <- jenks.tests(Jclassif) #requires the 'classInt' package
    df$class <- cut(df$sorted.values, unique(Jclassif$brks), labels=FALSE, include.lowest=TRUE) #the function unique() is used to remove non-unique breaks, should the latter be produced. This is done because the cut() function cannot break the values into classes if non-unique breaks are provided
    if(length(Jclassif$brks)!=length(unique(Jclassif$brks))){
      info <- ("The method has produced non-unique breaks, which have been removed. Please, check '...$classif$brks'")
    } else {info <- ("The method did not produce non-unique breaks.")}
    loop.res <- numeric(nrow(df))
    i <- 1
    repeat{
      i <- i+1
      loop.class <- classIntervals(df$sorted.values, i, style = "jenks")
      loop.test <- jenks.tests(loop.class)
      loop.res[i] <- loop.test[[2]]
      if(loop.res[i]>0.9999){
        break
      }
    }
    max.GoF.brks <- which.max(loop.res)
    plot(x=df$sorted.values, y=c(1:nrow(df)), type="b", main=paste0("Jenks natural breaks optimization; number of classes: ", n), sub=paste0("Goodness of Fit: ", round(test[[2]],4), ". Max GoF (", round(max(loop.res),4), ") with classes:", max.GoF.brks), ylim =c(0, nrow(df)+top.margin), cex=0.75, cex.main=0.95, cex.sub=0.7, ylab="observation index", xlab="value (increasing order)")
    abline(v=Jclassif$brks, lty=3, col="red")
    text(x=Jclassif$brks, y= max(nrow(df)) + dist, labels=sort(round(Jclassif$brks, 2)), cex=brks.cex, srt=90)
    results <- list("info"=info, "classif" = Jclassif, "breaks.max.GoF"=max.GoF.brks, "class.data" = df)
    return(results)
  }
  
  

# a function to retrieve taxonomy to accepted names and retrieve taxonomic hierarchy for a df with a column of taxonomic names
  
  gbif_tax <- function(df
                       , sppCol=1
                       , outFile="data/luGBIF.feather"
                       , kingType="Plantae"
                       , getCommon = FALSE
                       , targetRank = "Species"
                       ){
    
    tmpFile <- paste0(gsub(".feather","",outFile),"_temp.feather")
    
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
    
    assign("luRank",luRank,envir = .GlobalEnv)
    
    targetSort <- luRank %>%
      dplyr::filter(Rank == targetRank) %>%
      dplyr::pull(sort)
    
    alreadyDone01 <- if(file.exists(outFile)) read_feather(outFile) %>%
      dplyr::distinct(originalName) %>%
      dplyr::pull()
      
    alreadyDone02 <- if(file.exists(tmpFile)) read_feather(tmpFile) %>%
      dplyr::distinct(originalName) %>%
      dplyr::pull()
    
    alreadyDone <- c(get0("alreadyDone01"),get0("alreadyDone02"))
    
    toCheck <- df %>%
      dplyr::select(sppCol) %>%
      dplyr::distinct() %>%
      dplyr::pull()
    
    taxa <- tibble(originalName = setdiff(toCheck,alreadyDone)) %>%
      dplyr::filter(!grepl("BOLD:.*\\d{4}",originalName)
                    , !is.na(originalName)
                    ) %>%
      dplyr::mutate(taxa = gsub("dead|\\'|\\?| spp\\.| sp\\.|#|\\s^","",originalName)
                    , taxa = gsub(" x .*$| X .*$","",taxa)
                    , taxa = str_squish(taxa)
                    )
    
    if(length(taxa$taxa)>0){
      
      for (i in 1:length(taxa$taxa)){
        
        print(taxa$taxa[i])
        
        taxGBIF <- name_backbone(taxa$taxa[i], kingdom = kingType) %>%
          dplyr::mutate(originalName = taxa$originalName[i])
        
        taxGBIF <- if(sum(grepl("acceptedUsageKey",names(taxGBIF)))>0) {
          
          name_usage(taxGBIF$acceptedUsageKey,return="data")$data %>%
            dplyr::mutate(matchType = "Synonym") %>%
            dplyr::rename(usageKey = key
                          , status = taxonomicStatus
                          ) %>%
            dplyr::mutate(originalName = taxa$originalName[i])
          
        } else {
          
          taxGBIF
          
        }
          
        if(getCommon) taxGBIF$Common <- get_gbif_common(taxGBIF$usageKey)
        
        taxGBIF$Taxa <- taxGBIF %>%
          tidyr::pivot_longer(where(is.numeric),names_to = "key") %>%
          dplyr::mutate(key = map_chr(key,~gsub("Key","",.))
                        , key = str_to_sentence(key)
                        ) %>%
          dplyr::filter(key %in% luRank$Rank) %>%
          dplyr::left_join(luRank, by = c("key" = "Rank")) %>%
          dplyr::filter(sort <= targetSort) %>%
          dplyr::filter(sort == max(sort)) %>%
          dplyr::select(tolower(luRank$Rank[luRank$sort == .$sort])) %>%
          dplyr::pull()
        
        taxGBIF$Stamp <- Sys.time()
        
        if(file.exists(tmpFile)) {
          
          write_feather(taxGBIF %>%
                          dplyr::bind_rows(read_feather(tmpFile))
                        , paste0(gsub(".feather","",outFile),"_temp.feather")
                        )
          
        } else {
          
          write_feather(taxGBIF
                        , paste0(gsub(".feather","",outFile),"_temp.feather")
                        )
          
          }
        
      }
      
      # Clean up results
      read_feather(tmpFile) %>%
        {if(!file.exists(outFile)) (.) else (.) %>% dplyr::bind_rows(read_feather(outFile))} %>%
        dplyr::group_by(originalName) %>%
        dplyr::filter(Stamp == max(Stamp)) %>%
        dplyr::ungroup() %>%
        write_feather(outFile)
      
      file.remove(tmpFile)
      
    } else {
      
      {warning( "No taxa supplied" ) }
      
    }
    
  }
  
  # Find common name from GBIF (key is the gbif 'useagekey')
  get_gbif_common <- function(key) {
    
    print(key)
    
    commonNames <- name_usage(key)$data %>%
      dplyr::select(contains("Key")) %>%
      dplyr::select(where(is.numeric)) %>%
      tidyr::pivot_longer(1:ncol(.),names_to = "key") %>%
      dplyr::mutate(key = map_chr(key,~gsub("Key","",.))
                    , key = str_to_sentence(key)
                    ) %>%
      dplyr::filter(key %in% luRank$Rank) %>%
      dplyr::left_join(luRank, by = c("key" = "Rank")) %>%
      dplyr::filter(sort == max(sort)) %>%
      dplyr::pull(value) %>%
      name_usage(data="vernacularNames")
    
    df <- commonNames$data %>%
      dplyr::select(any_of(c("vernacularName","language","preferred")))
    
    hasAny <- nrow(df) > 0
    
    hasPreferred <- if("preferred" %in% names(df)) sum(df$preferred, na.rm = TRUE) > 0 else FALSE
    
    hasLanguage <- if("language" %in% names(df)) sum(df$preferred, na.rm = TRUE) > 0 else FALSE
    
    hasPreferredEng <- if(hasPreferred) df %>%
      dplyr::filter(preferred) %>%
      dplyr::filter(language == "eng") %>%
      nrow() %>%
      `>` (0) else FALSE
    
    hasEng <- if(hasLanguage) df %>%
      dplyr::filter(language == "eng") %>%
      nrow() %>%
      `>` (0) else FALSE
    
    if(hasPreferredEng) {
      
      df %>%
        dplyr::filter(preferred
                      , language == "eng"
                      ) %>%
        dplyr::pull(vernacularName) %>%
        paste0(collapse = ", ")
      
    } else if(hasEng) {
      
      df %>%
        dplyr::filter(language == "eng") %>%
        tidytext::unnest_tokens("common",vernacularName,token = "regex", pattern = ",|and",collapse = FALSE) %>%
        dplyr::mutate(common = gsub("^\\s|\\s$|etc","",common)) %>%
        dplyr::distinct(common) %>%
        dplyr::pull(common) %>%
        sort() %>%
        paste0(collapse = ", ")
      
    } else if(hasAny) {
      
      df %>%
        dplyr::count(language,vernacularName) %>%
        dplyr::arrange(desc(n)
                       , language
        ) %>%
        dplyr::slice(1) %>%
        dplyr::pull(vernacularName) %>%
        `[` (1)
      
    } else ""
    
  }
  
  # Add common name to existing taxonomic data frame
  add_gbif_common <- function(path = "data/luGBIF.feather") {
    
    gbifTaxDf <- read_feather(path) %>%
      #(if(testing) {. %>% dplyr::sample_n(5)} else {.}) %>%
      dplyr::mutate(Common = future_map_chr(key,get_gbif_common))
    
    write_feather(gbifTaxDf,path)
    
  }
  
  # Create polygon mask for filtering of patches
  make_aoi <- function(polygons
                       ,filterPolys = FALSE
                       ,filterPolysCol=NULL
                       ,polyBuffer
                       ,doMask = TRUE
                       ) {
    
    if(filterPolys != FALSE) polygons <- polygons %>%
        dplyr::filter(!!ensym(filterPolysCol) %in% filterPolys)
    
    polygons <- if(doMask) {
      
      polygons %>%
        dplyr::mutate(dissolve = 1) %>%
        dplyr::summarise(Include = n()) %>%
        st_cast() %>%
        st_buffer(polyBuffer) 
      
    } else {
      
      polygons %>%
        dplyr::mutate(dissolve = 1) %>%
        dplyr::summarise(Include = n()) %>%
        st_cast() %>%
        st_buffer(polyBuffer) %>%
        st_bbox() %>%
        st_as_sfc() %>%
        st_sf() %>%
        dplyr::mutate(Include = 1)
      
    }
    
    polygons <- polygons %>%
      st_transform(crs = 3577)
    
  }
  
  