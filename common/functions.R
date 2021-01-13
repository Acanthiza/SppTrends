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
                    , longText = paste0(var2, " occured ", text, " ", tolower(substr(colour,1,4)), " in ", var1, " than ",var1b)
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
  vec_to_sentence <- function(x,sep=",") {
    
    x[!is.na(x)] %>%
      paste(collapse = "JOINSRUS") %>%
      (function(x) if(sep == ";") {
        
        stringi::stri_replace_last_regex(x,"JOINSRUS", paste0(sep," and ")) %>%
          str_replace_all("JOINSRUS",paste0(sep," "))
        
      } else {
        
        stringi::stri_replace_last_regex(x,"JOINSRUS", " and ") %>%
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
    helper(x)
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
  