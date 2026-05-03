#' Fits a Random Forest of Interaction Trees
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
#' @description Estimates individualized treatment effects (ITEs) using Random Forest of Interaction Trees.
#' Works with randomized controlled trials (RCTs) and observational data.
#' Treatment variables may be binary, categorical, ordered, or continuous. For binary and survival outcomes, useRes = TRUE must be specified.
#' @param formula Formula to build CERFIT. Categorical predictors must be listed as a factor. e.g., Y ~ x1 + x2 | treatment.
#' For survival outcomes, the response may be specified as Surv(time, event) or Surv(start, stop, event).
#' @param data Data to grow a tree.
#' @param ntrees Number of trees to grow. This value should not be too small, as observations may not be well represented and averaging is insufficient, leading to unstable results. 
#' A value of 1000 is often recommended, but it can be reduced for smaller datasets or to speed up computation.
#' @param subset A logical vector that controls what observations are used to grow the forest.
#' The default value will use the entire data frame.
#' @param search Split search strategy. Options: “exhaustive” (evaluate all cut points) or “sss” (sigmoid approximation, experimental).
#' @param method Study type. "RCT" for randomized data, "observational" for observational data.
#' @param PropForm Method for estimating propensity scores (if method = "observational"). Options: "randomForest", "CBPS", "GBM", "HI". Not all options are compatible with all treatment types. See details.
#' @param split Impurity measure splitting statistic. Currently supports "t.test".
#' @param mtry Number of variables to consider at each split
#' @param nsplit Number of candidate cut points. If NULL, all possible cut points are considered. If an integer is provided, that many cut points are randomly selected from the set of possible cut points.
#' @param nsplit.random Default is FALSE. If TRUE, candidate cut points are chosen randomly.
#' @param minsplit Number of observations required to continue growing tree.
#' @param minbucket Number of observations required in each child node.
#' @param maxdepth Maximum depth of tree.
#' @param oob Logical, whether or not to use out-of-bag sample for predictions. Default is FALSE.
#' @param a Sigmoid approximation variable (for "sss" which is still under development).
#' @param sampleMethod Method to sample learning sample. Default is bootstrap. Subsample
#' takes a subsample of the original data. SubsamplebyID samples by an ID column and
#' uses all observations that have that ID. allData uses the entire data set
#' for every tree.
#' @param useRes Default is TRUE. If TRUE, fits the model using residuals from a regression of the response on covariates (excluding treatment).
#' Improves stability and accuracy. Linear regression is used for continuous responses, logistic regression for binary responses, Cox proportional hazards model are used for survival responses.
#' If response is binary or survival, useRes must be set equal to TRUE.
#' @param scale.y Logical, standardize y when creating splits (For "sss" to increase stability).
#' @param response Response type. Options are "auto", "continuous", "binary", and "survival".
#' When response = auto, the function automatically detects a survival response if the formula uses Surv(...), detects a binary response if the outcome is stored as a factor with exactly two levels, 
#' and otherwise treats the response as continuous. Thus, numeric continuous outcomes are recognized automatically, but binary outcomes coded numerically (for example, 0/1) are not automatically identified as binary
#' and should be specified with response = binary if that behavior is desired.
#' @param surv_resid Residual type used for survival responses. Options are "deviance" and "martingale".
#' @param surv_id Optional subject ID column name for survival models with time-dependent covariates in counting-process format, e.g. Surv(start, stop, event). Required for that format.
#' @return Returns a fitted CERFIT object which is a list containing the following elements:
#' \itemize{
#' \item randFor: The Random Forest of Interaction Trees.
#' \item trt.type: A string containing the treatment type of the data used to fit the model.
#' Can be binary, multiple, ordered, or continuous.
#' \item response.type: A string representing the response type of the data. Can be
#' binary or continuous.
#' \item useRes: A logical indicator that is TRUE if the model was fit on the
#' residuals of a linear model.
#' \item data: The data used to fit the model and includes the estimated propensity score if
#'  method was set to observational.
#' }
#' @details This function implements Random Forest of Interaction Trees proposed
#' in Su (2018), which is a modification of the Random Forest algorithm where
#' instead of a split being chosen to maximize prediction accuracy; each split
#' is chosen to maximize subgroup treatment heterogeneity. It chooses the best
#' split by maximizing the test statistic for \eqn{H_0: \beta_3=0} in the
#' following linear model:
#'
#' \eqn{Y_i = \beta_0 + \beta_1I(X_{ij} < c) + \beta_2I(Z = 1) + \beta_3I(X_{ij} < c)I(Z = 1) + \varepsilon_i}
#'
#' Where \eqn{X_{ij}} represents the splitting variable and Z = 1 represents
#' treatment. So, by maximizing the test statistic for \eqn{\beta_3} we are
#' maximizing the treatment difference between the nodes.
#'
#' The above equation only works when the data comes from a randomized controlled
#' trial, but we can modify it to give us unbiased estimates of treatment
#' effect in observational studies, as shown by Li et al. (2022). To do that we add propensity score into the
#' linear model.
#'
#'\eqn{Y_i = \beta_0 + \beta_1I(X_{ij} < c) + \beta_2I(Z = 1) + \beta_3I(X_{ij} < c)I(Z = 1) + \beta_4e_i + \varepsilon_i}
#'
#'Where \eqn{e_i} represents the propensity score. The CERFIT function will estimate
#'propensity score automatically when the method argument is set to observational.
#'
#'To control how this function estimates propensity score, you can use the
#'PropForm argument. Which can take four possible values randomForest, CBPS,
#' GBM and HI. randomForest uses the randomForest package to use a random forest
#' to estimate propensity score, CBPS uses Covariate Balancing Propensity Score
#' to estimate propensity score, GBM uses generalized boosted regression models
#' to estimate propensity score, and HI is the Hirano–Imbens generalized 
#' propensity score weighting for continuous treatments.
#' Some of these options only work for certain treatment types. 
#' See the full list below.
#' \itemize{
#'  \item binary: GBM, CBPS, randomForest
#'  \item categorical: GBM, CBPS
#'  \item ordered: GBM, CBPS
#'  \item continuous: CBPS, HI
#' }
#' Note: CBPS option supports up to four categories/levels of treatment, 
#' if there are more than four use GBM.
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
                    mtry=NULL, nsplit=NULL, nsplit.random=FALSE, minsplit=20, minbucket=round(minsplit/3), maxdepth=30,oob = FALSE,
                    a=50, sampleMethod=c('bootstrap','subsample','subsampleByID','allData'),
                    useRes=TRUE, scale.y=FALSE,
                    response=c("auto","continuous","binary","survival"),
                    surv_resid=c("deviance","martingale"),
                    surv_id=NULL)#
{
  sampleMethod <- match.arg(sampleMethod, c('bootstrap','subsample','subsampleByID','allData'))

  response <- match.arg(response, c("auto","continuous","binary","survival"))
  surv_resid <- match.arg(surv_resid, c("deviance","martingale"))
  if (!is.null(surv_id) && !is.character(surv_id)) {
    stop("surv_id must be a character string naming the subject ID column.", call. = FALSE)
  }

  # Detect survival response of the form Surv(time, event) or Surv(start, stop, event)
  lhs <- formula[[2]]
  is_surv <- is.call(lhs) && identical(as.character(lhs[[1]]), "Surv")
  surv_is_counting <- FALSE
  surv_time_vars <- NULL
  surv_event_var <- NULL

  # Build an internal response variable so downstream code can keep using all.vars(formula)[1]
  y_name <- all.vars(formula)[1]
  trt_name <- all.vars(formula)[length(all.vars(formula))]

  if (is_surv) {
    surv_args <- as.list(lhs)[-1]
    n_surv_args <- length(surv_args)

    if (n_surv_args == 2) {
      surv_time_vars <- as.character(surv_args[[1]])
      surv_event_var <- as.character(surv_args[[2]])
    } else if (n_surv_args == 3) {
      surv_time_vars <- c(as.character(surv_args[[1]]), as.character(surv_args[[2]]))
      surv_event_var <- as.character(surv_args[[3]])
      surv_is_counting <- TRUE
      if (is.null(surv_id)) {
        stop("surv_id must be provided when using time-dependent covariates with Surv(start, stop, event).", call. = FALSE)
      }
      if (!surv_id %in% names(data)) {
        stop("surv_id must name a column in data.", call. = FALSE)
      }
    } else {
      stop("Survival responses must be of the form Surv(time, event) or Surv(start, stop, event).", call. = FALSE)
    }

    # Covariates: everything except survival time variables, event indicator, and treatment
    covars <- setdiff(all.vars(formula)[-length(all.vars(formula))], c(surv_time_vars, surv_event_var))

    # Internal response variable name used for residuals
    y_name <- ".CERFIT_y"

    # Rebuild formula as: .CERFIT_y ~ x1 + ... + xp | treatment
    rhs <- if (length(covars) > 0) paste(covars, collapse = " + ") else "1"
    formula <- stats::as.formula(paste(y_name, "~", rhs, "|", trt_name))

    # If user didn't explicitly set response, treat Surv(...) as survival
    if (response == "auto") response <- "survival"
  }

  response.type <- "continuous"

  if (response == "survival") {
    response.type <- "survival"
  } else {
    response_vec <- data[[all.vars(formula)[1]]]
    if (is.factor(response_vec) && length(levels(response_vec)) == 2) {
      response.type <- "binary"
    }
  }

  response_print <- paste0(toupper(substring(response.type,first = 1,last = 1)),
                          substring(response.type,first = 2))
  cat(paste(response_print,"Response","\n"))

  if (useRes) {
    if (response.type == "survival") {
      # For survival: fit Cox model excluding treatment, then use Cox residuals as the working response
      covars <- all.vars(formula)[2:(length(all.vars(formula)) - 1)]
      rhs <- if (length(covars) > 0) paste(covars, collapse = " + ") else "1"

      if (surv_is_counting) {
        start_name <- surv_time_vars[1]
        stop_name <- surv_time_vars[2]
        event_name <- surv_event_var

        resformula <- stats::as.formula(
          paste0("survival::Surv(", start_name, ",", stop_name, ",", event_name, ") ~ ", rhs)
        )
        cox_fit <- survival::coxph(resformula, data = data, id = data[[surv_id]])

        # Preserve original survival outcome components
        data$yo_tstart <- data[[start_name]]
        data$yo_tstop <- data[[stop_name]]
        data$yo_delta <- data[[event_name]]
        data$yo_id <- data[[surv_id]]
      } else {
        time_name <- surv_time_vars[1]
        event_name <- surv_event_var

        resformula <- stats::as.formula(
          paste0("survival::Surv(", time_name, ",", event_name, ") ~ ", rhs)
        )
        cox_fit <- survival::coxph(resformula, data = data)

        # Preserve original survival outcome components
        data$yo_time <- data[[time_name]]
        data$yo_delta <- data[[event_name]]
      }

      # Residual choice: martingale or deviance
      if (surv_resid == "deviance") {
        eres <- stats::resid(cox_fit, type = "deviance")
      } else {
        eres <- stats::resid(cox_fit, type = "martingale")
      }

      # Store working response residuals in the internal response column
      data[[all.vars(formula)[1]]] <- as.numeric(eres)

    } else if (response.type == "binary") {
      resformula <- stats::as.formula(paste(all.vars(formula)[1], paste(all.vars(formula)[2:(length(all.vars(formula)) - 1)], collapse = " + "), sep = " ~ "))
      reslm <- stats::glm(resformula, data, family = stats::binomial)
      eres <- (as.numeric(data[[all.vars(formula)[1]]]) - 1) - stats::fitted(reslm)
      data$yo <- data[[all.vars(formula)[1]]]
      data[[all.vars(formula)[1]]] <- eres

    } else {
      resformula <- stats::as.formula(paste(all.vars(formula)[1], paste(all.vars(formula)[2:(length(all.vars(formula)) - 1)], collapse = " + "), sep = " ~ "))
      reslm <- stats::lm(resformula, data)
      eres <- stats::resid(reslm)
      data$yo <- data[[all.vars(formula)[1]]]
      data[[all.vars(formula)[1]]] <- eres
    }
  } else {
    if (response.type == "survival") {
      # Preserve original survival outcome components for reference
      if (surv_is_counting) {
        data$yo_tstart <- data[[surv_time_vars[1]]]
        data$yo_tstop <- data[[surv_time_vars[2]]]
        data$yo_delta <- data[[surv_event_var]]
        data$yo_id <- data[[surv_id]]
      } else {
        data$yo_time <- data[[surv_time_vars[1]]]
        data$yo_delta <- data[[surv_event_var]]
      }
      # Ensure internal response column exists (set to 0 if not residualizing)
      data[[all.vars(formula)[1]]] <- 0
    } else {
      data$yo <- data[[all.vars(formula)[1]]]
    }
  }
  trt_var <- all.vars(formula)[length(all.vars(formula))]
  TrT <- data[[trt_var]]
  trt_unique <- unique(TrT)
  trt.length <- length(trt_unique)
  if (trt.length < 2) stop("Only one treatment?", call. = FALSE)
  trt.type <- ifelse(trt.length == 2, "binary", "multiple")
  trt.type <- ifelse(is.ordered(TrT), "ordered", trt.type)
  trt.type <- ifelse(is.numeric(TrT) && trt.length > 10, "continuous", trt.type)
  trtlevels <- seq_len(trt.length)
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
        TrT_num <- as.numeric(TrT)
        propfun <- stats::lm(propformula, data = data[all.vars(formula)[-1]])
        prt <- as.numeric(stats::predict(propfun))
        sigm <- summary(propfun)$sigma
        prop <- stats::dnorm(TrT_num, mean = prt, sd = sigm)
        modhi <- stats::lm(TrT_num ~ 1)
        ps.num <- stats::dnorm(TrT_num, mean = as.numeric(modhi$fitted), sd = summary(modhi)$sigma)
        Iptw <- ps.num / prop
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
        prop <- matrix(NA,ncol= length(unique(TrT)),nrow=nrow(data))
        propfun <- twang::mnps(propformula,data=data[,all.vars(formula)[-1]],interaction.depth = 4,
                               stop.method = "es.max",estimand="ATE",verbose=FALSE,n.trees = 10000)
        pslist <- propfun$psList
        for(i in 1:length(unique(TrT))){
          prop[,i]<-unlist(pslist[[i]]$ps)
        }
        prop <- t(apply(prop, 1, cumsum))
        prop <- as.data.frame(prop)[,colSums(is.na(prop)) == 0]
        prop <- prop / prop[,length(unique(TrT))]
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
        prop <- matrix(NA,ncol= length(unique(TrT)),nrow=nrow(data))
        propfun <- twang::mnps(propformula,data=data[,all.vars(formula)[-1]],interaction.depth = 4,
                               stop.method = "es.max",estimand="ATE",verbose=FALSE,n.trees = 10000)
        pslist <- propfun$psList
        for(i in 1:length(unique(TrT))){
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
  if (response.type == "survival") {
    data$surv_resid_type <- surv_resid
    data$surv_is_counting <- surv_is_counting
    if (!is.null(surv_id)) data$surv_id_col <- surv_id
  }
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
    if(oob) {oob <- data[-obs.b,]} else {oob <- data[obs.b,]}
    tree.b <- growTree(formula=formula, data=sample.b,oob = oob, subset=subset, search=search, method=method, split=split,
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
