#' Calculate Variable Importance
#'
#' @param cerfit A fitted CERFIT object
#' @return Returns a named vector with the name of each predictor used to fit the CERFIT
#' object and its corresponding average minimal depth across all trees
#' @description Calculates the average minimal depth of each predictor used to fit
#' a CERFIT object. It calculates Variables importance by using a Variables average minimal depth.
#' variable's with a lower average minimal depth are more important.
#' @details  The depth of the root node is zero and if a variable does not appear
#' at any split in a tree it is assigned maxdepth + 1 for that tree.
#' @examples
#' fit <- CERFIT(Result_of_Treatment ~ sex + age + Number_of_Warts + Area + Time + Type | treatment,
#' data = warts,
#' ntrees = 30,
#' method = "RCT",
#' mtry = 2)
#' importance <- MinDepth(fit)
#' @export
MinDepth <- function(cerfit){  # need to given number of levels if observation
  cerfit <- cerfit$randFor
  Term<-cerfit[[1]]$tree$terms
  dataTemp<-all.vars(Term[[3]])
  vars<-dataTemp[-length(dataTemp)]

  mindepth <- rep(0, length(vars))
  for (t in seq_along(cerfit)) {
    intNodes <- nodeids(cerfit[[t]]$tree)[-nodeids(cerfit[[t]]$tree, terminal = TRUE)]
    varsInTree <- vars[unique(unlist(nodeapply(cerfit[[t]]$tree, ids=intNodes, FUN=function(n){split_node(n)$varid})))]
    varsAtNode <- unlist(nodeapply(cerfit[[t]]$tree, ids=intNodes, FUN=function(n){split_node(n)$varid}))
    #Root node should be counted as 0
    #depthAtNode <- table(unlist(lapply(intNodes, function(x) intersect(intNodes, nodeids(cerfit[[t]]$tree, from=x)))))-1
    #depthAtNode <- idDepth(cerfit[[t]]$tree)
    depthAtNode <- unlist(nodeapply(cerfit[[t]]$tree, ids = intNodes, info_node)) - 2
    treeDepth <- depth(cerfit[[t]]$tree)

    for (j in seq_along(vars)) {
      if (is.element(vars[j], varsInTree)) { #If variable is in tree
        mindepth[j]=mindepth[j]+min(depthAtNode[varsAtNode==j])
      } else {  #If variable not in tree, set mindepth to maximum depth+1
        mindepth[j]=mindepth[j]+treeDepth+1
      }
    }
  }
  mindepth <- mindepth/length(cerfit)
  names(mindepth) <- vars
  return(mindepth)
}
