#' Get predictions from a CERFIT object
#'
#' @param object An object of class CERFIT
#' @param data Sata used to build the tree
#' @param newdata New data to make predictions from
#' @param gridval For continuous treatment. Controls for what values of treatment to predict
#' @param prediction Return prediction using all trees ("overall") or using first i trees ("by iter")
#' @param type Choose what you want ro predict
#' @param alpha something
#' @param useRse THis will be removed
#' @param ... Additional Arguments
#' @export
predict.CERFIT <- function(object,data,newdata, gridval=NULL,
                           prediction=c("overall","by iter"),
                           type=c("response","ITE","node","opT"),
                           alpha=0.5,useRse=FALSE,...){

  #Return prediction using all trees ("overall") or using first i trees ("by iter")S
  prediction <- match.arg(prediction, c("overall","by iter"))

  type <- match.arg(type, c("response","ITE","node","opT"))
  cumMeanNA <- function(x){
    xTemp<-x;
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
  #ntrt <- length(utrt)
  # if grival is null, use the 10th quantile
  if(useRse){
    resformula <-  stats::as.formula(paste(all.vars(formulaTree)[1], paste(all.vars(formulaTree)[2:(length(all.vars(formulaTree))-1)], collapse=" + "), sep=" ~ "))
    reslm <- stats::lm(resformula,data)
    ylmp <- stats::predict(reslm,newdata)
    print("WHAT")
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
  } else{
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
    ntree<-length(x)
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
