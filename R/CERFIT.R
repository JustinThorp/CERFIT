#' Fits a Random Forest of Interactions Trees
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
#' @description Estimates an observations individualized treatment effect for RCT
#' and observational data. Treatment can be an binary, categorical, ordered, or continuous
#' variable. Currently if response is binary useRes must be set equal to TRUE.
#' @param formula Formula to build CERFIT.  Categorical predictors must be listed as a factor. e.g., Y ~ x1 + x2 | treatment
#' @param data Data to grow a tree.
#' @param ntrees Number of Trees to grow
#' @param subset A logical vector that controls what observations are used to grow the forest.
#' The default value will use the entire dataframe
#' @param search Method to search through candidate splits
#' @param method For observational study data, method="observational";for randomized study data, method="RCT".
#' @param PropForm Method to estimate propensity score
#' @param split Impurity measure splitting statistic
#' @param mtry Number of variables to consider at each split
#' @param nsplit Number of cut points selected
#' @param nsplit.random Logical: indicates if process to select cut points are random
#' @param minsplit Number of observations required to continue growing tree
#' @param minbucket Number of observations required in each child node
#' @param maxdepth Maximum depth of tree
#' @param a Sigmoid approximation variable (for "sss" which is still under development)
#' @param sampleMethod Method to sample learning sample. Default is bootstrap. Subsample
#' takes a subsample of the original data. SubsamplebyID samples by an ID column and
#' uses all observations that have that ID. allData uses the entire data set
#' for every tree.
#' @param useRes Logical indicator if you want to fit the CERFIT model to
#' the residuals from a linear model
#' @param scale.y Logical, standardize y when creating splits (For "sss" to increase stability)
#' @return Returns a fitted CERFIT object which is a list with the following elements
#' \itemize{
#' \item RandFor: The Random forest of interaction trees
#' \item trt.type: A string containing the treatment type of the data used to fit the model.
#' Cant be binary, multiple, ordered or continuous.
#' \item response.type: A string representing the response type of the data. Can be
#' binary or continuous.
#' \item useRes: A logical indicator that is TRUE if the model was fit on the
#' residuals of a linear model
#' \item data: The data used to fit the model also contains the propensity score if
#'  method was set to observational}
#' @details This function implements Random Forest of Interaction Trees proposed
#' in Su (2018). Which is a modification of the Random Forest algorithm where
#' instead of a split being chosen to maximize prediction accuracy each split
#' is chosen to maximized subgroup treatment heterogeneity. It chooses the best
#' split by maximizing the test statistic for \eqn{H_0: \beta_3=0} in the
#' following linear model
#'
#' \eqn{Y_i = \beta_0 + \beta_1I(X_{ij} < c) + \beta_2I(Z = 1) + \beta_3I(X_{ij} < c)I(Z = 1) + \varepsilon_i}
#'
#' Where \eqn{X_{ij}} represents the splitting variable and Z = 1 represents
#' treatment. So, by maximizing the  test statistic for \eqn{\beta_3} we are
#' maximizing the treatment difference between the nodes.
#'
#' The above equation only works when the data comes from a randomized controlled
#' trial. But we can modify it to gives us unbiased estimates of treatment
#' effect in observational studies Li et al. (2022). To do that we add propensity score into the
#' linear model.
#'
#'\eqn{Y_i = \beta_0 + \beta_1I(X_{ij} < c) + \beta_2I(Z = 1) + \beta_3I(X_{ij} < c)I(Z = 1) + \beta_4e_i + \varepsilon_i}
#'
#'Where \eqn{e_i} represents the propensity score. The CERIT function will estimate
#'propensity score automatically when the method argument is set to observational.
#'
#'To control how this function estimates propensity score you can use the
#'PropForm argument. Which can take four possible values randomForest, CBPS,
#' GBM and HI. randomForest uses the randomForest package to use a random forest
#' to estimate propensity score, CBPS uses Covariate balancing propensity score
#' to estimate propensity score GBM uses generalized boosted regression models
#' to estimate propensity score, and HI is for continuous treatment and
#' estimates the general propensity score. Some of these options only work
#' for certain treatment types. Full list below
#' \itemize{
#'  \item binary: GBM, CBPS, randomForest
#'  \item categorical: GBM, CBPS
#'  \item ordered: GBM, CBPS
#'  \item continuous: CBPS, HI
#' }
#'
#' @references
#' \itemize{
#' \item Li, Luo, et al. Causal Effect Random Forest of
#' Interaction Trees for Learning Individualized Treatment Regimes with
#'  Multiple Treatments in Observational Studies. Stat, 2022,
#'  https://doi.org/10.1002/sta4.457.
#' \item Su, X., Peña, A., Liu, L., & Levine, R. (2018). Random forests of interaction trees for estimating individualized treatment effects in randomized trials.
#' Statistics in Medicine, 37(17), 2547- 2560.
#' \item G. W. Imbens, The role of the propensity score in estimating dose-response
#' functions., Biometrika, 87 (2000), pp. 706–710.
#' \item G. Ridgeway, D. McCarey, and A. Morral, The twang package: Toolkit for
#' weighting and analysis of nonequivalent groups, (2006).
#' \item A. Liaw and M. Wiener, Classification and regression by randomforest, R
#' News, 2 (2002), pp. 18–22}
#' @examples
#' fit <- CERFIT(Result_of_Treatment ~ sex + age + Number_of_Warts + Area + Time + Type | treatment,
#' data = warts,
#' ntrees = 30,
#' method = "RCT",
#' mtry = 2)
#'
#' @export
### Grows a random forest ###
# Res is for fitting the residuals
CERFIT <- function( formula, data, ntrees, subset = NULL,search=c("exhaustive","sss"),
                    method=c("RCT","observational"), PropForm=c("randomForest","CBPS","GBM", "HI"),
                    split=c("t.test"),
                    mtry=NULL, nsplit=NULL, nsplit.random=FALSE, minsplit=20, minbucket=round(minsplit/3), maxdepth=30,
                    a=50, sampleMethod=c('bootstrap','subsample','subsampleByID','allData'),
                    useRes=TRUE, scale.y=FALSE)#
{
  sampleMethod <- match.arg(sampleMethod, c('bootstrap','subsample','subsampleByID','allData'))
  if (missing(formula)) stop("A formula must be supplied.", call. = FALSE)
  if (missing(data)) data <- NULL
  response <- data[[all.vars(formula)[1]]]
  response.type = "continous"
  if(is.factor(response) & length(levels(response)) == 2){
    response.type = "binary"
    #useRes = FALSE # Residual dont work with binary response right now
  }
  response_print <- paste0(toupper(substring(response.type,first = 1,last = 1)),
                          substring(response.type,first = 2))
  cat(paste(response_print,"Response","\n"))
  if(useRes){
    if (response.type == "binary") {
      resformula<- stats::as.formula(paste(all.vars(formula)[1], paste(all.vars(formula)[2:(length(all.vars(formula))-1)], collapse=" + "), sep=" ~ "))
      reslm <- stats::glm(resformula,data,family = stats::binomial)
      eres <- (as.numeric(data[[all.vars(formula)[1]]]) - 1) - stats::fitted(reslm)
      data$yo <- data[[all.vars(formula)[1]]]
      data[[all.vars(formula)[1]]] <- eres
    } else {
      resformula<- stats::as.formula(paste(all.vars(formula)[1], paste(all.vars(formula)[2:(length(all.vars(formula))-1)], collapse=" + "), sep=" ~ "))
      reslm <- stats::lm(resformula,data)
      eres <- stats::resid(reslm)
      data$yo <- data[[all.vars(formula)[1]]]
      data[[all.vars(formula)[1]]] <- eres
    }
  } else {
    data$yo <- data[[all.vars(formula)[1]]]
  }
  TrT <- data[all.vars(formula)[length(all.vars(formula))]]
  trt.length<-nrow(unique(TrT))
  if (trt.length<2) stop("Only one treatment?", call. = FALSE)
  trt.type <- ifelse(trt.length==2,"binary","multiple")
  trt.type <- ifelse(is.ordered(TrT[[1]]),"ordered",trt.type)
  trt.type <- ifelse(trt.length>10, "continuous", trt.type)
  trtlevels<-c(1:trt.length)
  trttype_print <- paste0(toupper(substring(trt.type,first = 1,last = 1)),
                            substring(trt.type,first = 2))
  cat("Treatment Levels: ")
  cat(paste0(trtlevels),"\n")
  cat(paste(trttype_print,"Treatment","\n"))


  if(method=="observational"){
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
      if(PropForm == "GBM") {
        prop <- matrix(NA,ncol= length(unique(TrT[[1]])),nrow=nrow(data))
        propfun <- twang::mnps(propformula,data=data[,all.vars(formula)[-1]],interaction.depth = 4,
                               stop.method = "es.max",estimand="ATE",verbose=FALSE,n.trees = 10000)
        pslist <- propfun$psList
        for(i in 1:length(unique(TrT[[1]]))){
          prop[,i]<-unlist(pslist[[i]]$ps)
        }
        prop <- t(apply(prop, 1, cumsum))
        prop <- as.data.frame(prop)[,colSums(is.na(prop)) == 0]
        prop <- prop / prop[,length(unique(TrT[[1]]))]
        #names(prop) <- TrT_splits[!is.na(TrT_splits)]
        Iptw <- twang::get.weights(propfun,stop.method = "es.max",estimand="ATE")
        #Iptw<- rep(1,nrow(data))
        #Iptw <- sum(TrT)/length(TrT)*TrT/prop+sum(1-TrT)/length(TrT)*(1-TrT)/(1-prop)
        #Iptw <-TrT/prop+(1-TrT)/(1-prop)
        Iptw <- truncquant(Iptw,q=0.9)
      } else if (PropForm == "CBPS") {
        #data[,all.vars(formula)[length(all.vars(formula))]]<-as.factor(data[,all.vars(formula)[length(all.vars(formula))]])
        propfun <- CBPS::CBPS(propformula, data = data[,all.vars(formula)[-1]],ATT=FALSE,method = "exact")#
        prop <- propfun$fitted.values
        print(prop)
        Iptw <- propfun$weights
        levels(data[,all.vars(formula)[length(all.vars(formula))]])<-c(1:trt.length)
      } else if (PropForm == "old") {
        prop <- matrix(NA,ncol= length(unique(TrT[[1]])),nrow=nrow(data))
        propfun <- twang::mnps(propformula,data=data[,all.vars(formula)[-1]],interaction.depth = 4,
                               stop.method = "es.max",estimand="ATE",verbose=FALSE,n.trees = 10000)
        pslist <- propfun$psList
        for(i in 1:length(unique(TrT[[1]]))){
          prop[,i]<-unlist(pslist[[i]]$ps)
        }
        prop <- t(apply(prop, 1, cumsum))
        prop <- as.data.frame(prop)[,colSums(is.na(prop)) == 0]
        #names(prop) <- TrT_splits[!is.na(TrT_splits)]
        Iptw <- twang::get.weights(propfun,stop.method = "es.max",estimand="ATE")
        #Iptw<- rep(1,nrow(data))
        #Iptw <- sum(TrT)/length(TrT)*TrT/prop+sum(1-TrT)/length(TrT)*(1-TrT)/(1-prop)
        #Iptw <-TrT/prop+(1-TrT)/(1-prop)
        Iptw <- truncquant(Iptw,q=0.9)
      }
    }  else stop("Please specify a propensity score method: randomForest or CBPS or GBM", call. = FALSE)
  } else if (method=="RCT") {
    prop <- rep(1,nrow(data))#rep("none",nrow(data)) # for observational no prop need
    Iptw<- rep(1,nrow(data))}
  if(!exists("Iptw")) {
    stop("Not able to estimate Propenisty Score. \n Check that your function arguments are correct")
  }
  #data[,all.vars(formula)[length(all.vars(formula))]]<-as.numeric(as.character(data[,all.vars(formula)[length(all.vars(formula))]]))
  data[,all.vars(formula)[-length(all.vars(formula))]] <- sapply(data[,all.vars(formula)[-length(all.vars(formula))]],as.numeric)
  data$iptw <- Iptw
  data$prop <- prop
  #data <- cbind(data,prop)
  #return(data)
  #Construct random forest
  randFor <- lapply(1:ntrees,function(b){
    if(b%%10==0){cat(paste0("Tree Number: ",b,"\n"))}
    #print(paste0("Tree Number: ",b))
    obs.b <- switch(sampleMethod,
                    bootstrap = sample.int(nrow(data), size=nrow(data), replace=TRUE, prob=data$iptw), #inverse weighting in boostrapping
                    subsample = sample.int(nrow(data), size=round(nrow(data)*0.632), replace=FALSE,prob=data$iptw), # stratified sampling
                    #subsampleByID = {nIds <- length(unique(data[[idVar]]))
                    #unlist(lapply(sample(unique(data[[idVar]]), size=round(nIds*0.632), replace=FALSE),
                    #              function(x){which(data[[idVar]] == x)}))},
                    allData = 1:nrow(data))
    sample.b <- data[obs.b,]
    tree.b <- growTree(formula=formula, data=sample.b, subset=subset, search=search, method=method, split=split,
                       mtry=mtry, nsplit=nsplit, nsplit.random=nsplit.random, minsplit=minsplit, minbucket=minbucket, maxdepth=maxdepth, a=a,
                       scale.y=scale.y, useRes=useRes, trtlevels=trtlevels,response.type = response.type)#, useRpart=useRpart, minpvalue=minpvalue, corstr=corstr)
    list(tree=tree.b,cases=sort(unique(obs.b)))
  })
  trt <- data[[all.vars(formula)[length(all.vars(formula))]]]
  #print(length(trt))
  #print(nrow(data))
  #print(length(all.vars(formula)))
  data[[all.vars(formula)[length(all.vars(formula))]]] <- as.numeric(as.character(trt))
  object <- list(randFor = randFor,trt.type = trt.type,
                 response.type = response.type,
                 useRes = useRes,
                 data = data)
  class(object) <- "CERFIT"
  return(object)
}
# Having issues in mutiple treatment where in partition only a single treatment
# is present in the data
