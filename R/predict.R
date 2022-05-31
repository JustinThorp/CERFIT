#' Get predictions from a CERFIT object
#'
#' @param object A fitted CERFIT object
#' @param newdata New data to make predictions from. IF not provided will make predictions
#' on training data
#' @param gridval For continuous treatment. Controls for what values of treatment to predict
#' @param prediction Return prediction using all trees ("overall") or using first i trees ("by iter")
#' @param type Choose what value you wish to predict. Response will predict the response.
#' ITE will predict the Individualized treatment effect. Node will predict the node. And opT
#' will predict the optimal treatment for each observation.
#' @param alpha For continuous treatment it is the mixing parameter for the elastic
#' net regularization in each node. When equal to 0 it is ridge regression and
#' when equal to 1 it is lasso regression.
#' @param ... Additional Arguments
#' @return The return value depends of the type argument. If type is response the function
#' will return a matrix with n rows and the number of columns equal to the level of treatment.
#' If type is ITE then it returns a matrix with n rows and a number of columns equal to
#' one minus the levels of treatment. And if type is opT then it returns a matrix with n
#' rows and two columns. With the first column denoting the optimal treatment and
#' the second column denoting the optimal response.
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
  response.type <- object$response.type
  treatment.type <- object$trt.type
  object <- object$randFor
  type <- match.arg(type, c("response","ITE","node","opT"))
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
  ## should add warnings here if gridbalue beyond min or max utrt
  ntrt <- length(utrt)
  # if grival is null, use the 10th quantile
  if(useRse == TRUE & response.type == "continous"){
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
  if(length(utrt)<=20){ ## if less than 20 unique treatments/levels using unique treatments
    ntrt=length(utrt)
    gridval<-utrt
  } else if(is.null(gridval)) { # if more than 20, and gridval is null, use percentiles at 5% increment
    gridval <- stats::quantile(utrt, prob = seq(0, 1, length = 21))
    ntrt<-length(gridval)-1
  } else {
    ntrt<-length(gridval)}
  print(gridval)

  if(type!="opT"){
    predictMat <- lapply(lapply(object, "[[" , "tree"), predictTree, newdata=newdata,gridval=gridval,ntrt=ntrt,type=type,LB=LB,UB=UB,alpha=alpha)
    ypre<- do.call(cbind,predictMat)
    #yp<- lapply(1:ntrt,function(i,k) k[,seq(i, by = ntrt, length = NCOL(ypre) / ntrt)],k=ypre)
    ypre<- lapply(1:ntrt,function(i,k) k[,seq(i, NCOL(ypre), by = ntrt)], k=ypre)
    y.pre<- t(matrix(unlist(lapply(ypre,rowMeans,na.rm=TRUE)), ncol=NROW(newdata),byrow = TRUE))
    y.pre<-y.pre+ylmp
    #y.pre: by row observation, each column is the corresponding predition for 1 treatment.
  } else if (type == "opT" && treatment.type != "continous"){
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
      yname<-paste("y",utrt,"-y",utrt[1],sep="")
      ite<- y.pre-y.pre[,1]
      colnames(ite) <- c(yname)
      return(ite[,-1])
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
      Ypre[[i]]<-t(apply(ypre[[i]],1,cumMeanNA))
    }
    cumypre<-t(matrix(unlist(Ypre),ncol=NROW(newdata),byrow = TRUE))
    ntree<-length(object)
    cumypre.l<- lapply(seq(1,(ntrt*ntree),by=ntree),function(i,k) k[,i:(i+ntree-1)], k=cumypre)
    print(cumypre.l)
    if(type=="response"){
      yname<-paste("ycum",utrt,sep="")
      names(cumypre.l) <- yname
      return(cumypre.l)}
    if(type=="ITE")  {
      cumite<-as.list(NA)
      for(i in 1:ntrt){
        cumite[[i]]<- cumypre.l[[i]]-cumypre.l[[1]]
      }
      yname<-paste("ycum",utrt,"-ycum",utrt[1],sep="")
      names(cumite)<-yname
      print(yname)
      return(cumite[[-1]])
    }}

}
