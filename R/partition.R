partition<- function(vars, y, trt, propensity, subset, search, method, split, nsplit, nsplit.random,
                     minsplit, minbucket, a, scale.y, useSearch, useOptim,trtlevels,response.type){#, allVars
  if (sum(subset) < minsplit) {return(NULL)}
  vars <- vars[subset,,drop=FALSE]
  y <- y[subset]
  trt<- trt[subset]
  if (length(unique(trt)) < 2) {return(NULL)}
  
  
  # if(length(trtlevels) > 2 & length(trtlevels) < 10 & method != "RCT") { # 10/7/25 JM fix typo & !is.ordered(trt)) {
  #   propensity <- propensity[subset,]
  # } else {
  #   propensity <- propensity[subset]
  # }
  
  # --- Normalize shape of 'propensity' before subsetting (handles GBM vector/data.frame/matrix cases)
  if (is.data.frame(propensity)) {
    # If GBM output was wrapped in a 1-column data.frame, flatten it
    if (ncol(propensity) == 1L) {
      propensity <- as.numeric(propensity[[1L]])
    } else {
      # Multi-arm GBM (multiple treatment levels) → convert to matrix
      propensity <- as.matrix(propensity)
    }
  } else if (is.matrix(propensity) && ncol(propensity) == 1L) {
    # 1-column matrix → flatten to numeric vector
    propensity <- as.numeric(propensity[, 1L])
  }
  
  # --- Subset correctly depending on treatment structure
  if (length(trtlevels) > 2 & length(trtlevels) < 10 & method != "RCT") {
    # multi-arm observational → subset rows + all columns
    propensity <- propensity[subset, , drop = FALSE]
  } else {
    # binary or RCT → subset rows only
    propensity <- propensity[subset]
  }
  
  
  trt.length <- length(trtlevels)
  trt_num_unique <- length(unique(trt))
  is_continuous_trt <- is.numeric(trt) && !is.factor(trt) && !is.ordered(trt) && trt_num_unique > 10
  
  if (method != "RCT") {
    if (is.ordered(trt) || is_continuous_trt) {
      # Treat ordered and numeric-continuous treatments the same way during split evaluation
      trt_num <- as.numeric(trt)
      trt_vals <- sort(unique(trt_num))
      if (length(trt_vals) < 2) {return(NULL)}
      ran <- sample(seq_len(length(trt_vals) - 1L), 1L)
      cutoff_trt <- trt_vals[ran]
      trt <- ifelse(trt_num <= cutoff_trt, 1, 0)

      # For ordered multi-level GPS matrices, align propensity with the temporary dichotomization
      if (is.matrix(propensity) && ncol(propensity) > 1L) {
        propensity <- rowSums(propensity[, seq_len(ran), drop = FALSE])
      }

    } else if (trt.length > 2 & trt.length < 10) {
      ## if less than 10 treatments/levels
      ran <- sample(unique(trt), 2)
      keep <- trt == ran[1] | trt == ran[2]
      vars <- vars[keep, , drop = FALSE]
      y <- y[keep]
      if (is.matrix(propensity)) {
        propensity <- propensity[keep, , drop = FALSE]
      } else {
        propensity <- propensity[keep]
      }
      trt <- trt[keep]
      trt <- ifelse(trt == ran[1], 1, 0)

      # If multi-arm and matrix, grab matching column
      if (is.matrix(propensity)) {
        col_id <- match(as.character(ran[1]), colnames(propensity))
        if (is.na(col_id)) {
          suppressWarnings(col_id <- as.integer(ran[1]))
        }
        if (is.na(col_id) || col_id < 1L || col_id > ncol(propensity)) {
          stop("Could not align propensity column with sampled treatment level")
        }
        propensity <- propensity[, col_id]
      }
    }
  }
  else {
    if (is.ordered(trt) || is_continuous_trt) {
      # Chooses the split point for ordered or numeric-continuous treatment
      trt_num <- as.numeric(trt)
      trt_vals <- sort(unique(trt_num))
      if (length(trt_vals) < 2) {return(NULL)}
      ran <- sample(seq_len(length(trt_vals) - 1L), 1L)
      cutoff_trt <- trt_vals[ran]
      trt <- ifelse(trt_num <= cutoff_trt, 1, 0)
    } else if (trt.length > 2 & trt.length < 10) {
      ## if less than 10 treatments/levels
      ran <- sample(unique(trt), 2)
      keep <- trt == ran[1] | trt == ran[2]
      vars <- vars[keep, , drop = FALSE]
      y <- y[keep]
      trt <- trt[keep]
      trt <- ifelse(trt == ran[1], 1, 0)
    }
  }


  if (NROW(vars) < 2*minbucket) {return(NULL)}
  if (length(unique(y))==1) {return(NULL)}
  stats<- cutoff<- breakLeft<-NA
  findStats<-sapply(vars,function(x){
    x_is_factor <- is.factor(x)
    x_levels <- if (x_is_factor) levels(x) else NULL
    if (search=="exhaustive" && !is.null(nsplit) && nsplit.random) {
      xTemp <- ordinalize(x, y, sortCat=FALSE)
    } else {
      xTemp <- ordinalize(x, y, sortCat=TRUE)
      }
    x <- xTemp$x
    #If all x values the same, do not check optimal split
    if (abs(max(x) - min(x)) > 1e-8) {
      #The SSS partition deals with problems when there is a very small number of observations
      #Use exhaustive search in this case (or set minsplit >= 5)
      if (search=="sss") { #leave sss here for now
        print("sss not ready")
      } else if (search=="exhaustive") { #current codes only work exhaustive search
        cutpts <- findCutpts(x, minbucket)
        #z <- matrix(x,ncol = length(x))[rep(1, length(cutpts)), ] < cutpts
        if (is.null(nsplit)) {
          nsplit<- length(cutpts)
        }
        #Take nsplit cutpoints (if applicable)
        if (!is.null(nsplit) && !is.null(cutpts) && length(cutpts) > 1) {
          #If nsplit.random is TRUE, take nsplit cutpts randomly.  Otherwise, take nsplit cutpts equally spread out across cutpts
          if (!nsplit.random & length(cutpts) > nsplit) { #if not random select nsplit cut
            cutpts <- unique(cutpts[seq(1, length(cutpts), length.out=nsplit)])
          } else {
            cutpts <- sort(sample(cutpts, min(c(nsplit, length(cutpts))), replace=FALSE))
          }
        }
        #
        #It is possible (unlikely) no cutpoint can satisfy minbucket
        if (!is.null(cutpts)) {
          #print(list(y,x,trt,cutpts,method,propensity,minbucket,response.type))
          mod <- find_split(y=y, x=x, trt=trt,cutpts=cutpts,
                            method=method,propensity = propensity,
                            minbucket=minbucket,response_type = response.type)
          if (!is.na(mod$stat)) {
            stats <- mod$stat
            if (x_is_factor) {
              cutoff <- "factor"
              breakLeft <- rep(NA_integer_, length(x_levels))
              names(breakLeft) <- x_levels
              left_lvls <- colnames(xTemp$cutToLvl)[xTemp$cutToLvl <= mod$cutoff]
              right_lvls <- colnames(xTemp$cutToLvl)[xTemp$cutToLvl > mod$cutoff]
              breakLeft[x_levels %in% left_lvls] <- 1L
              breakLeft[x_levels %in% right_lvls] <- 2L
              if (all(is.na(breakLeft)) || length(unique(stats::na.omit(breakLeft))) <= 1) {
                stop("Did not find correct cutpoints")
              }
            }
            else {cutoff <- mod$cutoff; breakLeft<-NA}
          }
        }

      }
      else {stop("Unexpected search")}
    }
    return(c(stats,cutoff,breakLeft))})
  #If randomly picking a subset of categories, do not sort by mean.  Would be more likely to select variables when sorted
  #return(findStats)
  #If each candidate variable cannot be split (e.g. cannot satisfy minbucket), return null
  if (all(is.na(findStats[1,]))) {return(NULL)}
  best_idx <- which.max(findStats[1,])
  if (identical(findStats[2, best_idx], "factor")) {
    # Index is used for categorical variable splits
    return(partysplit(varid = as.integer(colnames(findStats)[best_idx]),
                      index = as.integer(findStats[3, best_idx]),
                      info = list(stats = findStats[1, ])))
  } else {

    #Breaks is used for continuous variable splits
    #print(as.integer(colnames(findStats)[which.max(findStats[1,])]))
    return(partysplit(varid = as.integer(colnames(findStats)[best_idx]),
                      breaks = findStats[2, best_idx],
                      info = list(stats = findStats[1, ])))
  }
}
