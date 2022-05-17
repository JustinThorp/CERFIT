partition<- function(vars, y, trt, propensity, subset, search, method, split, nsplit, nsplit.random,
                     minsplit, minbucket, a, scale.y, useSearch, useOptim,trtlevels,response.type){#, allVars
  if (sum(subset) < minsplit) {return(NULL)}
  vars <- vars[subset,,drop=FALSE]
  y <- y[subset]
  trt<- trt[subset]
  if (length(unique(trt)) < 2) {return(NULL)}
  if(length(trtlevels) > 2 & length(trtlevels<10) & method != "RCT") { #& !is.ordered(trt)) {
    propensity <- propensity[subset,]
  } else {
    propensity <- propensity[subset]
  }
  trt.length<-length(trtlevels)
  if (method != "RCT") {
    if (is.ordered(trt)) {
      # Chosses the split point for ordered treatment
      ran <- sample(1:(length(propensity) - 2),1) # fix this so it chooses right split point
      # and stuff
      propensity <- propensity[,ran]
      trt <- ifelse(trt <= ran,1,0)
    } else if (trt.length>2 & trt.length < 10) {
      ## if less than 10 treatments/levels
      ran<- sample(unique(trt),2)
      vars<-subset(vars,trt==ran[1] | trt==ran[2])
      y<-subset(y,trt==ran[1] | trt==ran[2])
      propensity<-subset(propensity,trt==ran[1] | trt==ran[2])
      trt<-subset(trt,trt==ran[1] | trt==ran[2])
      trt<-ifelse(trt==ran[1],1,0)
      propensity<-propensity[,ran[1]] # need to make sure trt is levels as propensity nameorders
    }
  } else {
    if (is.ordered(trt)) {
      # Chooses the split point for ordered treatment
      trt <- as.numeric(trt)
      ran <- sample(min(trt):(max(trt) - 2),1) # fix this so it chooses right split point
      # and stuff
      trt <- ifelse(trt <= ran,1,0)
    } else if (trt.length>2 & trt.length < 10) {
      ## if less than 10 treatments/levels
      ran<- sample(unique(trt),2)
      vars<-subset(vars,trt==ran[1] | trt==ran[2])
      y<-subset(y,trt==ran[1] | trt==ran[2])
      #propensity<-subset(propensity,trt==ran[1] | trt==ran[2])
      trt<-subset(trt,trt==ran[1] | trt==ran[2])
      trt<-ifelse(trt==ran[1],1,0)
      #propensity<-propensity[,ran[1]]
    }
  }


  if (NROW(vars) < 2*minbucket) {return(NULL)}
  if (length(unique(y))==1) {return(NULL)}
  stats<- cutoff<- breakLeft<-NA
  findStats<-sapply(vars,function(x){
    x_factor_check <- as.numeric(x)
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
            if (is.factor(x_factor_check)) {
              cutoff<- "factor"
              breakLeft <- rep(NA, length(levels(x)))
              breakLeft[levels(x) %in% colnames(xTemp$cutToLvl)[xTemp$cutToLvl <= mod$cutoff]]=1L
              breakLeft[levels(x) %in% colnames(xTemp$cutToLvl)[xTemp$cutToLvl > mod$cutoff]]=2L
              print(breakLeft)
              if (all(is.na(breakLeft)) & length(unique(breakLeft))<=1) {stop("Did not find correct cutpoints")}
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
  if (inherits(findStats[2,which.max(findStats[1,])],"factor")) {
    #Index is used for categorical variable splits
    print("factor")
    return(partysplit(varid=as.integer(colnames(findStats)[which.max(findStats[1,])]),
                      index=findStats[3,which.max(findStats[1,])],
                      info=list(stats=findStats[1,])))
  } else {

    #Breaks is used for continuous variable splits
    #print(as.integer(colnames(findStats)[which.max(findStats[1,])]))
    return(partysplit(varid=as.integer(colnames(findStats)[which.max(findStats[1,])]),
                      breaks=findStats[2, which.max(findStats[1,])],
                      info=list(stats=findStats[1,])))
  }
}

