### Convert factors to numerical value. ###
ordinalize <- function(x, y, sortCat=TRUE){
  if (is.factor(x)) {
    x <- factor(x) #Remove factors not listed
    #One can randomly assign a category a distinct numerical value
    if (!sortCat) {
      cutToLvl <- t(sample.int(length(levels(x))))
      colnames(cutToLvl)=levels(x)
    } else {
      #For binary, sort data by proportion in class 1.  For continuous, sort by means
      if (is.factor(y)) {
        cutToLvl <- prop.table(table(y,x),2)[1,,drop=FALSE]
      } else {cutToLvl <- t(vapply(levels(x), function(z){mean(y[x==z])}, numeric(1)))}
    }
    #Convert lvls to numerical value. Slow method. Make this faster later.
    xTemp <- rep(NA,length(x))
    for (lvls in levels(x)) {
      xTemp[x==lvls] <- cutToLvl[colnames(cutToLvl)==lvls]
    }
  } else {
    xTemp <- x
    cutToLvl <- NULL
  }
  return(list(x=xTemp, cutToLvl=cutToLvl))
}
