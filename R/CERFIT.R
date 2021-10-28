# Functions ---------------

### Load Functions ###
#rm(list=ls(all=TRUE))
#library(partykit)
#library(parallel)
#library(pROC)
#library(CBPS)
#library(randomForest)
#library(twang)
#library(glmnet)
#library(sandwich)
#library(Rcpp)
#library(RcppArmadillo)
##sourceCpp("\\src\\find_split.cpp")
## usethis namespace: start
#' @useDynLib CERFIT, .registration = TRUE
## usethis namespace: end
## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
#' @importFrom partykit partysplit
#' @importFrom partykit party
#' @importFrom partykit partynode
#' @importFrom partykit kidids_split
#' @importFrom partykit nodeids
#' @importFrom partykit fitted_node
#' @importFrom partykit as.constparty
#' @importFrom partykit nodeapply
#' @importFrom partykit split_node
#' @importFrom partykit info_node
#' @importFrom grid depth
#' @importFrom stats complete.cases
#' @importFrom stats terms


#'
#' @param formula Formula to build CERFIT.  Categorical predictors must be listed as a factor. e.g., Y ~ x1 + x2 | treatment
#' @param data Data to grwo a tree.
#' @param search Method to search throught candidate splits
#' @param method For observational stuy data, method="observation";for randomized study data, method="RCT".
#' @return The fitted CERFIT model.
#' @examples
#' add(1, 1)
#' add(10, 1)
#' @export
### Grows a random forest ###
# Res is for fitting the residuals
CERFIT <- function( formula, data, ntrees, subset=NULL, search=c("exhaustive","sss"),
                    method=c("RCT","observation"), PropForm=c("randomForest","CBPS","GBM", "HI"),
                    split=c("t.test"),
                    mtry=NULL, nsplit=NULL, nsplit.random=TRUE, minsplit=20, minbucket=round(minsplit/3), maxdepth=30,
                    a=50, sampleMethod=c('bootstrap','subsample','subsampleByID','learning'),
                    useRes=TRUE, scale.y=FALSE)#
{
  sampleMethod <- match.arg(sampleMethod, c('bootstrap','subsample','subsampleByID','learning'))
  if (missing(formula)) stop("A formula must be supplied.", call. = FALSE)
  if (missing(data)) data <- NULL

  response <- data[[all.vars(formula)[1]]]
  response.type = "continous"
  if(is.factor(response) & length(levels(response)) == 2){
    response.type = "binary"
    #useRes = FALSE # Residual dont work with binary response right now
  }
  print(paste(response.type,"response"))
  if(useRes){
    resformula<- stats::as.formula(paste(all.vars(formula)[1], paste(all.vars(formula)[2:(length(all.vars(formula))-1)], collapse=" + "), sep=" ~ "))
    reslm <- stats::lm(resformula,data)
    eres <- stats::resid(reslm)
    data$yo <- data[[all.vars(formula)[1]]]
    data[[all.vars(formula)[1]]] <- eres
  }
  TrT <- data[all.vars(formula)[length(all.vars(formula))]]
  trt.length<-nrow(unique(TrT))
  if (trt.length<2) stop("Only one treatment?", call. = FALSE)
  trt.type <- ifelse(trt.length==2,"binary","multiple")
  trt.type <- ifelse(is.ordered(TrT[[1]]),"ordered",trt.type)
  #trt.type <- ifelse(trt.length>10, "continuous", trt.type)
  trtlevels<-c(1:trt.length)
  print(trtlevels)
  print(paste(trt.type,"Treatment"))


  if(method=="observation"){
    propformula <- stats::as.formula(paste(all.vars(formula)[length(all.vars(formula))], paste(all.vars(formula)[2:(length(all.vars(formula))-1)], collapse=" + "), sep=" ~ "))
    if(trt.type=="continuous"){
      if(PropForm=="CBPS"){
        propfun <- CBPS::CBPS(propformula, data = data[,all.vars(formula)[-1]],ATT=FALSE,method = "exact")#
        prop <- propfun$fitted.values
        Iptw <- propfun$weights
      } else if(PropForm=="HI") {
        propfun <- stats::lm(propformula,data=data[all.vars(formula)[-1]])
        prt <- stats::predict(propfun)
        sigm <- summary(propfun)$sigma
        prop <- stats::dnorm(TrT,prt,sigm)
        modhi = stats::lm(TrT~1)
        ps.num = stats::dnorm((TrT-modhi$fitted)/(summary(modhi))$sigma,0,1)
        Iptw=ps.num/prop
      }
    } else if(trt.type=="binary") {
      if (PropForm=="GBM") {
        propfun <- twang::ps(propformula,data=data[,all.vars(formula)[-1]],interaction.depth = 4, stop.method = "es.max",estimand="ATE",verbose=FALSE,n.trees = 10000)
        prop <- propfun$ps
        Iptw<- twang::get.weights(propfun,stop.method = "es.max",estimand="ATE")
      } else if (PropForm=="CBPS") {
        propfun <- CBPS::CBPS(propformula, data = data[,all.vars(formula)[-1]],ATT=FALSE,method = "exact")#
        prop <- propfun$fitted.values
        Iptw <- propfun$weights
      } else if (PropForm=="randomForest") {
        propfun<- suppressWarnings(randomForest::randomForest(propformula,data=data[all.vars(formula)[-1]]))
        prop <- propfun$predicted
        Iptw <- sum(TrT)/length(TrT)*TrT/prop+sum(1-TrT)/length(TrT)*(1-TrT)/(1-prop)
        #Iptw <-TrT/prop+(1-TrT)/(1-prop)
        Iptw <- truncquant(Iptw[[1]],q=0.9)
      }
    } else if (trt.type=="multiple"){
      if(PropForm=="GBM") {
        data[,all.vars(formula)[length(all.vars(formula))]]<-as.factor(data[,all.vars(formula)[length(all.vars(formula))]])
        propfun <- twang::mnps(propformula,data=data[,all.vars(formula)[-1]],interaction.depth = 4, stop.method = "es.max",estimand="ATE",verbose=FALSE,n.trees = 10000)
        pslist<-propfun$psList
        prop<-matrix(NA,ncol=trt.length,nrow=nrow(data))
        for(i in 1:trt.length){
          prop[,i]<-unlist(pslist[[i]]$ps)
        }
        colnames(prop)<-levels(data[,all.vars(formula)[length(all.vars(formula))]])
        levels(data[,all.vars(formula)[length(all.vars(formula))]])<-c(1:trt.length)
        Iptw <- twang::get.weights(propfun,stop.method = "es.max",estimand="ATE")
      } else if (PropForm=="CBPS" & trt.length<5 ) {
        data[,all.vars(formula)[length(all.vars(formula))]]<-as.factor(data[,all.vars(formula)[length(all.vars(formula))]])
        propfun <- CBPS::CBPS(propformula, data = data[,all.vars(formula)[-1]],ATT=FALSE,method = "exact")#
        prop <- propfun$fitted.values
        Iptw <- propfun$weights
        levels(data[,all.vars(formula)[length(all.vars(formula))]])<-c(1:trt.length)

      }
    } else if(trt.type == "ordered") {
      prop <- matrix(NA,ncol= length(unique(TrT[[1]])),nrow=nrow(data))
      order_length <- seq(0,length(unique(TrT[[1]])))
      i <- 1
      TrT_splits <- rep(NA,length(order_length))
      for (TrT_split in order_length){
        TrT_temp <- TrT[[1]] <= TrT_split
        if (length(unique(TrT_temp)) == 1) next else TrT_splits[i] <- TrT_split
        propformula_temp <- stats::as.formula(paste("TrT_temp",
                         paste(all.vars(formula)[2:(length(all.vars(formula))-1)],
                               collapse = "+"),
                         sep = "~"))
        propfun <- CBPS::CBPS(propformula_temp,ATT=FALSE,method = "exact",
                        data = data[,all.vars(formula)[c(-1,-length(all.vars(formula)))]])
        #propfun<- suppressWarnings(randomForest(x = data[all.vars(formula)[-length(all.vars(formula))][-1]],
        #                                        y = as.factor(TrT_temp)))
        prop[,i] <- propfun$fitted.values
        i <- i + 1
      }
      #data[,all.vars(formula)[length(all.vars(formula))]] <- as.factor(as.numeric(data[,all.vars(formula)[length(all.vars(formula))]]))
      prop <- as.data.frame(prop)[,colSums(is.na(prop)) == 0]
      names(prop) <- TrT_splits[!is.na(TrT_splits)]
      #data[,all.vars(formula)[length(all.vars(formula))]]<-as.factor(as.numeric(data[,all.vars(formula)[length(all.vars(formula))]]))
      print(class(data[,all.vars(formula)[length(all.vars(formula))]]))
      propfun <- twang::mnps(propformula,data=data[,all.vars(formula)[-1]],interaction.depth = 4,
                    stop.method = "es.max",estimand="ATE",verbose=FALSE,n.trees = 10000)
      Iptw <- twang::get.weights(propfun,stop.method = "es.max",estimand="ATE")
      #Iptw<- rep(1,nrow(data))
      #Iptw <- sum(TrT)/length(TrT)*TrT/prop+sum(1-TrT)/length(TrT)*(1-TrT)/(1-prop)
      #Iptw <-TrT/prop+(1-TrT)/(1-prop)
      Iptw <- truncquant(Iptw[[1]],q=0.9)
    }  else stop("Please specify a propensity score method: randomForest or CBPS or GBM", call. = FALSE)
  } else if (method=="RCT") {
    prop <- rep(1,nrow(data))#rep("none",nrow(data)) # for observational no prop need
    Iptw<- rep(1,nrow(data))}
  #data[,all.vars(formula)[length(all.vars(formula))]]<-as.numeric(as.character(data[,all.vars(formula)[length(all.vars(formula))]]))
  data[,all.vars(formula)[-length(all.vars(formula))]] <- sapply(data[,all.vars(formula)[-length(all.vars(formula))]],as.numeric)
  data$iptw <- Iptw
  data$prop <- prop
  #data <- cbind(data,prop)
  #return(data)
  #Construct random forest
  randFor <- lapply(1:ntrees,function(b){
    if(b%%10==0){print(paste0("Tree Number: ",b))}
    #print(paste0("Tree Number: ",b))
    obs.b <- switch(sampleMethod,
                    bootstrap = sample.int(nrow(data), size=nrow(data), replace=TRUE, prob=data$iptw), #inverse weighting in boostrapping
                    subsample = sample.int(nrow(data), size=round(nrow(data)*0.632), replace=FALSE,prob=data$iptw), # stratified sampling
                    #subsampleByID = {nIds <- length(unique(data[[idVar]]))
                    #unlist(lapply(sample(unique(data[[idVar]]), size=round(nIds*0.632), replace=FALSE),
                    #              function(x){which(data[[idVar]] == x)}))},
                    learning = 1:nrow(data))
    sample.b <- data[obs.b,]
    tree.b <- growTree(formula=formula, data=sample.b, subset=subset, search=search, method=method, split=split,
                       mtry=mtry, nsplit=nsplit, nsplit.random=nsplit.random, minsplit=minsplit, minbucket=minbucket, maxdepth=maxdepth, a=a,
                       scale.y=scale.y, useRes=useRes, trtlevels=trtlevels,response.type = response.type)#, useRpart=useRpart, minpvalue=minpvalue, corstr=corstr)
    list(tree=tree.b,cases=sort(unique(obs.b)),trt.type = trt.type)
  })
  class(randFor) <- "CERFIT"
  return(randFor)
}
# Having issues in mutiple treatment where in partition only a single treatment
# is present in the data
