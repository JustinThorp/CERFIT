### Return predicted values from a tree. ###
predictTree <- function(tree, newdata=tree$data, gridval, LB, UB, ntrt,type="response",alpha ){
  parse_tree_formula <- function(form) {
    form_chr <- paste(deparse(form), collapse = " ")
    form_chr <- gsub("\\s+", " ", form_chr)

    split_tilde <- strsplit(form_chr, "~", fixed = TRUE)[[1]]
    if (length(split_tilde) != 2) {
      stop("Unable to parse fitted CERFIT tree formula.", call. = FALSE)
    }

    rhs_full <- trimws(split_tilde[2])
    rhs_parts <- strsplit(rhs_full, "|", fixed = TRUE)[[1]]
    if (length(rhs_parts) != 2) {
      stop("Unable to identify predictors from fitted CERFIT tree formula.", call. = FALSE)
    }

    x_part <- trimws(rhs_parts[1])
    trt_part <- trimws(rhs_parts[2])

    predictors <- trimws(unlist(strsplit(x_part, "\\+")))
    predictors <- predictors[nzchar(predictors)]

    list(
      predictors = predictors,
      treatment = trt_part
    )
  }
  da <- stats::fitted(tree)
  colnames(da)[2:5]<-c("y","Trt","prop","ww")
  ufit<-sort(unique(da[["(fitted)"]]))

  formulaTree <- stats::formula(tree$terms)
  formula_info <- parse_tree_formula(formulaTree)
  predictors <- formula_info$predictors
  treatment <- formula_info$treatment

  required_cols <- unique(c(predictors, treatment))
  if (!all(required_cols %in% colnames(newdata))) {
    missing_cols <- setdiff(required_cols, colnames(newdata))
    stop(
      paste(
        "newdata is missing columns required for tree prediction:",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  newdata_tree <- newdata[, required_cols, drop = FALSE]

  tree_tmp <- tree
  tree_tmp$terms <- stats::terms(
    stats::as.formula(
      paste("~", paste(required_cols, collapse = " + "))
    )
  )
  nodesNewdata <- stats::predict(tree_tmp, newdata=newdata_tree, type="node")
  if (length(nodesNewdata) != nrow(newdata_tree)) {
    stop("Tree node assignment failed for one or more rows in newdata.", call. = FALSE)
  }
  Y.min<-ifelse(min(da[,2])<0,2*min(da[,2]),min(da[,2])/2)
  Y.max<-ifelse(max(da[,2])<0,max(da[,2])/2,2*max(da[,2]))
  #if (ntrt<=10){ #for binary and multiple trt, ignore gridval
  if(ntrt<=10){
    pred<- lapply(split(da ,list(da[["(fitted)"]],da[,3])), function(da){
      ytemp<-try(stats::weighted.mean(da[,2],da[,length(da)],na.rm=T))#This is only using the first propensity (fixed)
      if(inherits(ytemp,"try-error")) {
        return(NA)
        } else {
          return(ytemp)
          }
    })
    nodepred<- cbind(ufit,t(matrix(unlist(pred), ncol = length(ufit), byrow = TRUE)))
  } else {
    pred<- lapply(split(da ,da[["(fitted)"]]), function(da){
      Trt<-da$Trt
      x<-cbind(Trt,Trt^2,Trt^3)
      lambdas <- 10^seq(5, -3, by = -.1)
      fit <- try(glmnet::cv.glmnet(x, da$y,  family = "gaussian", alpha = alpha, lambda = lambdas,nfolds =10),silent=TRUE)
      if (inherits(fit, "try-error")){
        fit2<-try(glmnet::cv.glmnet(x, jitter(da$y), family = "gaussian", alpha = alpha, lambda = lambdas, nfolds =10),silent=TRUE)
        if (inherits(fit2, "try-error")) {
          return(NA)} else {fit<-fit2}
      }

      Trt<-gridval
      ext<-Trt>max(da[,3])|Trt<min(da[,3])
      nd<-cbind(gridval,gridval^2,gridval^3)
      ytemp <- stats::predict(fit, newx = nd, s=fit$lambda.min)
      ytemp[!ext]=ifelse(ytemp[!ext]>Y.max,Y.max,ytemp[!ext])##avoid extrem values
      ytemp[!ext]=ifelse(ytemp[!ext]<Y.min,Y.min,ytemp[!ext])#mean(da[,2])
      ytemp[ext]=ifelse(ytemp[ext]>Y.max,NA,ytemp[ext])
      ytemp[ext]=ifelse(ytemp[ext]<Y.min,NA,ytemp[ext])
      if (type!="opT") {
        return(ytemp)
      }else {
        #top<-gridval[which.max(ytemp)]
        #yop<-max(ytemp)
        lengthout<-5
        B <- seq(LB, UB, length.out=lengthout)
        opY<-Y.min; opTrt <- NA
        pref<-function(Trt){
          trtt<-cbind(Trt,Trt^2,Trt^3)
          yp<- stats::predict(fit, newx = trtt, s=fit$lambda.min)
          return(yp)}
        for (b in 1:(lengthout-1)) {
          fit.tmp <-  suppressWarnings(stats::optimize(pref, lower=B[b], upper=B[b+1], maximum=TRUE))
          if (is.na(fit.tmp$objective)) {
            opY<-opY
            opTrt<-opTrt
          } else {
            if (!is.nan(fit.tmp$objective) && fit.tmp$objective > opY && fit.tmp$objective < Y.max ) {
              opY <- fit.tmp$objective
              opTrt <- fit.tmp$maximum
            }
          } }
        return(cbind(opTrt,opY))}
    })
    nodepred<- cbind(ufit,matrix(unlist(pred), ncol = length(pred[[1]]), byrow = TRUE))
  }
  if(type=="opT" && ntrt  > 10) {
    ntrt<-2
    }
  predictions<-as.data.frame(cbind(nodesNewdata,matrix(NA,ncol=ntrt,nrow=nrow(newdata))))
  matched_rows <- match(predictions$nodesNewdata,nodepred[,1])
  if (nrow(nodepred) == 0 || all(is.na(matched_rows))) {
    stop("No terminal-node predictions could be matched to newdata in predictTree().", call. = FALSE)
  }
  predictions[,2:(ntrt+1)] <- nodepred[matched_rows,2:(ntrt+1), drop = FALSE]
  return(predictions[,2:(ntrt+1)])
}
