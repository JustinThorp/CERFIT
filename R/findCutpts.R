###Finds candidate cutpoints ###
findCutpts <- function(x, minbucket) {
  nX <- length(x)
  x_sort <- sort(x,method = "quick")
  #Due to ties, it's possible minbucket cannot be satisfied
  if (x_sort[minbucket]==rev(x_sort)[minbucket]) {
    cutpts=NULL}
  else {
    #Cutpoints considered must satisfy minbucket
    cutpts <- unique(x_sort[minbucket:(nX-minbucket+1)])

    if(length(cutpts)==1){stop(paste0("Only 1 cutpt??? ", cutpts, x))}
    cutpts <- (cutpts[1:(length(cutpts)-1)]+cutpts[2:length(cutpts)])/2
  }
  return(cutpts)
}
