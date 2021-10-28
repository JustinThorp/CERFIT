truncquant <- function(prop,q=0.9){
  qmax <- stats::quantile(prop,q)
  qmin <- stats::quantile(prop,1-q)
  prop[prop >= qmax] <- qmax
  prop[prop<=qmin] <- qmin
  return(prop)
}
