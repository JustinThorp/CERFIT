#' Get predictions from a CERFIT object
#'
#' @param object A fitted CERFIT object
#' @param newdata New data to make predictions from. If not provided will make predictions
#' on training data.
#' @param gridval For continuous treatments with more than 15 unique values, the
#' treatment values at which to predict. Must contain more than 10 values;
#' shorter grids return an error. If \code{NULL}, the 0th, 5th, ..., 100th
#' percentiles of the unique observed treatment values (21 values) are used. When the
#' treatment has 15 or fewer unique values, \code{gridval} is ignored and
#' predictions are made at each observed treatment value.
#' @param prediction Return prediction using all trees ("overall") or using first i trees ("by iter").
#' @param type Choose which value you wish to predict: 'response' will predict the potential outcome.
#' 'ITE' will predict the individualized treatment effect. And 'opT'
#' will predict the optimal treatment for each observation.
#' @param alpha For continuous treatment. It is the mixing parameter for the elastic
#' net regularization in each node. When equal to 0 it is ridge regression and
#' when equal to 1 it is lasso regression.
#' @param ... Additional Arguments
#' @return The return value depends on the \code{prediction} and \code{type} arguments.
#'
#' With \code{prediction = "overall"}: if \code{type = "response"}, a matrix with
#' n rows and one column per treatment level. If \code{type = "ITE"}, a matrix
#' with n rows and one column per contrast against the reference (lowest)
#' treatment level, i.e. the number of levels minus one. If \code{type = "opT"},
#' a matrix with n rows and two columns: the optimal treatment and the optimal
#' response.
#'
#' With \code{prediction = "by iter"}: a named list of matrices, each with n rows
#' and one column per tree. Column \code{j} is the estimate averaged over the
#' first \code{j} trees, skipping trees that give no prediction for that
#' observation, so the last column equals the \code{"overall"} result.
#' \itemize{
#' \item \code{type = "response"}: one matrix per treatment level (per
#' \code{gridval} value for a continuous treatment), named \code{"ycum<level>"}.
#' \item \code{type = "ITE"}: one matrix per contrast against the reference
#' level, named \code{"ycum<level>-ycum<reference>"}. A list is returned even
#' for a binary treatment, where it has length one.
#' }
#' \code{type = "opT"} is not supported with \code{prediction = "by iter"}.
#' @section Reproducibility with continuous treatments:
#' For a continuous treatment (a numeric treatment with more than 10 distinct
#' values), each terminal node fits a cubic dose-response curve with
#' \code{glmnet::cv.glmnet()}, which assigns cross-validation folds at random.
#' As a result, repeated calls to \code{predict()} on the same fitted object can
#' return slightly different values, and separate calls with
#' \code{type = "response"}, \code{"ITE"}, or \code{"opT"} will not exactly
#' reconcile with one another. To make predictions reproducible, call
#' \code{set.seed()} immediately before each \code{predict()} call, using the
#' same seed for calls whose results will be compared.
#' Alternatively, predict the responses once and compute the ITEs directly as
#' \code{resp[, -1] - resp[, 1]}, which is identical to \code{type = "ITE"}.
#'
#' Binary, multi-arm categorical, and ordered treatments, including numeric
#' treatments with 10 or fewer distinct values, do not use this path; their
#' predictions are deterministic, so a seed is needed only when fitting the
#' model with \code{CERFIT()}.
#' @examples
#' fit <- CERFIT(Result_of_Treatment ~ sex + age + Number_of_Warts + Area + Time + Type | treatment,
#' data = warts,
#' ntrees = 30,
#' method = "RCT",
#' mtry = 2)
#' ite <- predict(fit,type = "ITE")
#' @export
predict.CERFIT <- function(object,newdata = NULL, gridval=NULL,
                           prediction=c("overall","by iter"),
                           type=c("response","ITE","node","opT"),
                           alpha=0.5,...){

  #Return prediction using all trees ("overall") or using first i trees ("by iter")S
  prediction <- match.arg(prediction, c("overall","by iter"))
  useRse <- object$useRes
  data <- object$data
  if (is.null(newdata)) newdata <- object$data

  # Coerce newdata columns to match the stored training-data classes
  common_vars <- intersect(names(newdata), names(data))
  for (v in common_vars) {
    train_col <- data[[v]]
    new_col <- newdata[[v]]

    if ((is.integer(train_col) || is.numeric(train_col)) && is.factor(new_col)) {
      tmp <- suppressWarnings(as.numeric(as.character(new_col)))
      if (!any(is.na(tmp) & !is.na(new_col))) {
        if (is.integer(train_col)) {
          newdata[[v]] <- as.integer(tmp)
        } else {
          newdata[[v]] <- tmp
        }
      }
    } else if ((is.factor(train_col) || is.ordered(train_col)) && !is.factor(new_col)) {
      newdata[[v]] <- factor(new_col, levels = levels(train_col), ordered = is.ordered(train_col))
    } else if (is.character(train_col) && !is.character(new_col)) {
      newdata[[v]] <- as.character(new_col)
    } else if (is.logical(train_col) && !is.logical(new_col)) {
      newdata[[v]] <- as.logical(new_col)
    } else if (is.numeric(train_col) && !is.numeric(new_col) && !is.factor(new_col)) {
      newdata[[v]] <- as.numeric(new_col)
    } else if (is.integer(train_col) && !is.integer(new_col) && !is.factor(new_col)) {
      newdata[[v]] <- as.integer(new_col)
    }
  }

  # Note: "node" is accepted for forward compatibility but is not yet implemented.
  # It is excluded from the package documentation. If developed in the future,
  # a corresponding return block should be added in the prediction == "overall"
  # and prediction == "by iter" sections below.

  response.type <- object$response.type
  treatment.type <- object$trt.type
  object <- object$randFor
  type <- match.arg(type, c("response","ITE","node","opT"))
  # JDM 10/18/2026 - Add stop when prediction == "by iter" && type == "opT"; not supported
  if (prediction == "by iter" && type == "opT") {
    stop('type = "opT" is not supported with prediction = "by iter".', call. = FALSE)
  }
  cumMeanNA <- function(x){
    xTemp <- x;
    xTemp[is.na(xTemp)] <- 0
    cumsum(xTemp)/cumsum(!is.na(x))
    }
  #utrt<- sort(unique(c(fitted(x[[1]]$tree)[,3],fitted(x[[2]]$tree)[,3],fitted(x[[3]]$tree)[,3])))
  formulaTree <- stats::formula(object[[1]]$tree$terms)
  treatment <- all.vars(formulaTree)[length(all.vars(formulaTree))]
  utrt<-sort(unique(data[[treatment]]))
  LB<-min(data[[treatment]])
  UB<-max(data[[treatment]])
  qu<-seq(LB,UB,length.out = 6)
  ## add a statement warning if gridvalue beyond the LB and UB
  ## should add warnings here if gridvalue beyond min or max utrt
  ntrt <- length(utrt)
  # if grival is null, use the 10th quantile
  if(useRse == TRUE & response.type == "continuous"){
    resformula <-  stats::as.formula(paste("yo", paste(all.vars(formulaTree)[2:(length(all.vars(formulaTree))-1)], collapse=" + "), sep=" ~ "))
    reslm <- stats::lm(resformula,data)
    ylmp <- stats::predict(reslm,newdata)
    #print("WHAT")
  } else if(useRse == TRUE & response.type == "binary") {
    resformula <-  stats::as.formula(paste("yo", paste(all.vars(formulaTree)[2:(length(all.vars(formulaTree))-1)], collapse=" + "), sep=" ~ "))
    reslm <- stats::glm(resformula,data,family = stats::binomial())
    ylmp <- stats::predict(reslm,newdata,type = "response")
    print(length(ylmp))
  } else {
    ylmp<-rep(0,nrow(newdata))
  }
  if(length(utrt)<=15){ ## if less than or equal to 15 unique treatments/levels using unique treatments
    ntrt=length(utrt)
    gridval<-utrt
  } else if(is.null(gridval)) { # if more than 15, and gridval is null, use percentiles at 5% increment
    gridval <- stats::quantile(utrt, prob = seq(0, 1, length = 21))
    #ntrt<-length(gridval)-1 # JDM 10/17/2026 - This was causing an error when gridval was NULL and the treatment type was continuous. It was setting ntrt to 20 instead of 21, which caused an error in the predictTree function. I changed it to length(gridval) to fix the error.
    ntrt <- length(gridval)
  } else {
    ntrt<-length(gridval)}
  # JDM 10/17/2026
  # predictTree() switches to the categorical method when ntrt <= 10,
  # which would ignore gridval for a continuous treatment
  if (treatment.type == "continuous" && ntrt <= 10) {
    stop("For continuous treatments, gridval must contain more than 10 values.",
         call. = FALSE)
  }
  print(gridval)

  if(type!="opT"){
    predictMat <- lapply(lapply(object, "[[" , "tree"), predictTree, newdata=newdata,gridval=gridval,ntrt=ntrt,type=type,LB=LB,UB=UB,alpha=alpha)
    ypre<- do.call(cbind,predictMat)
    #yp<- lapply(1:ntrt,function(i,k) k[,seq(i, by = ntrt, length = NCOL(ypre) / ntrt)],k=ypre)
    ypre<- lapply(1:ntrt,function(i,k) k[,seq(i, NCOL(ypre), by = ntrt)], k=ypre)
    y.pre<- t(matrix(unlist(lapply(ypre,rowMeans,na.rm=TRUE)), ncol=NROW(newdata),byrow = TRUE))
    y.pre<-y.pre+ylmp
    #y.pre: by row observation, each column is the corresponding predition for 1 treatment.
  } else if (type == "opT" && treatment.type != "continuous"){
    predictMat<-lapply(lapply(object , "[[" , "tree"), predictTree, newdata=newdata,gridval=gridval,ntrt=ntrt,type="opT",  LB=LB,UB=UB,alpha=alpha)
    #ntrt<-2
    ypre<- do.call(cbind,predictMat)
    ypre<- lapply(1:ntrt,function(i,k) k[,seq(i, NCOL(ypre), by = ntrt)], k=ypre)
    y.pre<- t(matrix(unlist(lapply(ypre,rowMeans,na.rm=TRUE)), ncol=NROW(newdata),byrow = TRUE))
    y.pre<- y.pre + ylmp
    t.opt <- max.col(y.pre)
    y.opt <- apply(y.pre, 1, max, na.rm = TRUE)
    #topt<-as.matrix(ypre[[1]])
    #yopt<-as.matrix(ypre[[2]])
    #y.opt<-rowMeans(yopt,na.rm = T)+ylmp
    #t.opt<-rowMeans(topt,na.rm = T)
    y.pre<- cbind(t.opt,y.opt)
  } else {
    predictMat<-lapply(lapply(object , "[[" , "tree"), predictTree, newdata=newdata,gridval=gridval,ntrt=ntrt,type="opT",  LB=LB,UB=UB,alpha=alpha)
    ntrt<-2
    ypre<- do.call(cbind,predictMat)
    ypre<- lapply(1:ntrt,function(i,k) k[,seq(i, NCOL(ypre), by = ntrt)], k=ypre)
    topt<-as.matrix(ypre[[1]])
    yopt<-as.matrix(ypre[[2]])
    y.opt<-rowMeans(yopt)+ylmp
    t.opt<-rowMeans(topt)
    y.pre<- cbind(t.opt,y.opt)
  }

  yname<-NA
  if (prediction=="overall") {
    if(type=="response") {
      resp <- y.pre
      yname<- paste("y=",gridval,sep="")
      colnames(resp) <- yname
      return(resp)}
    if(type=="ITE") { #using the first level or smallest value as reference group
      yname <- paste("y", gridval[-1], "-y", gridval[1], sep = "")
      ite <- y.pre - y.pre[,1]
      ite <- ite[,-1, drop = FALSE]
      colnames(ite) <- yname
      return(ite)
    }

    if(type=="opT") {
      yname<-c("opTreat","opResponse")
      opTY<-y.pre
      colnames(opTY) <- c(yname)
      return(opTY)
    }
  }
  else if(prediction=="by iter"){
    Ypre<-as.list(NA)
    for(i in 1: ntrt){
      #Ypre[[i]]<-t(apply(ypre[[i]],1,cumMeanNA)) # JDM 10/17/2026 - response omits the ylmp offset
      Ypre[[i]] <- t(apply(ypre[[i]], 1, cumMeanNA)) + ylmp
    }
    cumypre<-t(matrix(unlist(Ypre),ncol=NROW(newdata),byrow = TRUE))
    ntree<-length(object)
    cumypre.l<- lapply(seq(1,(ntrt*ntree),by=ntree),function(i,k) k[,i:(i+ntree-1)], k=cumypre)
    print(cumypre.l)
    if(type=="response"){
      #yname<-paste("ycum",utrt,sep="") JDM 10/17/2026 - wrong labels for continuous treatments
      yname <- paste("ycum", gridval, sep = "")
      names(cumypre.l) <- yname
      return(cumypre.l)}
    if(type=="ITE")  {
      cumite<-as.list(NA)
      for(i in 1:ntrt){
        cumite[[i]]<- cumypre.l[[i]]-cumypre.l[[1]]
      }
      #yname<-paste("ycum",utrt,"-ycum",utrt[1],sep="")
      yname <- paste("ycum", gridval, "-ycum", gridval[1], sep = "")
      names(cumite)<-yname
      print(yname)
      #return(cumite[[-1]]) JDM 10/17/2026 - only will work with two arms, return maxtrix
      return(cumite[-1]) # will return a list one matrix per contrast against the reference arm
    }}

}
