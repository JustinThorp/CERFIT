growTree <- function(formula, data, subset=NULL, search=c("exhaustive","sss"),
                     method=c("RCT","observational"),
                     split=c("t.test", "pvalue"),#, "gini", "entropy", "information"),
                     mtry=NULL, nsplit=NULL, nsplit.random=TRUE, minsplit=20, minbucket=round(minsplit/3), maxdepth=20,
                     a=50, useRes, scale.y=FALSE, trtlevels,response.type)
{
  search <- match.arg(search,c("exhaustive","sss"))
  method <- match.arg(method,c("RCT","observational"))
  split <- match.arg(split,c("t.test", "pvalue"))
  stopifnot(is.logical(nsplit.random), is.logical(scale.y), is.logical(useRes))#, is.logical(useRpart))
  if (is.numeric(nsplit) && !nsplit.random && nsplit < 5) {"Selecting <5 ordered splits may yield unexpected results"}
  response <- all.vars(formula)[1]
  if(grepl("\\|", as.character(formula)[3])){
    treatment <- trimws(strsplit(as.character(formula)[3], "\\|")[[1]][2], which="both")
  } else {stop("Please specify the treatment in formula")}

  Propensity <- "prop"
  Iptw <- "iptw"
  data <- data[c(all.vars(formula)[-1], response, Iptw,Propensity)] #Rearrange data so that response comes last

  if (!all(complete.cases(data[-length(data)])) & !is.null(subset)) { paste0("Specifying subset with missing data can yield unexpected results") }
  data <- data[complete.cases(data[-length(data)]),]
  # The reason for not checking the last column for NA's is that the last column can be a
  # data.frame(?) and would through an error. Since the last column is propensity any
  # NA's are caused by NA is other columns so this should be fine


  if (is.null(mtry)){mtry <- length(all.vars(formula[[3]]))-1} #defult mtry=p

  #if(is.factor(data[[response]])){data[[response]]=as.numeric(data[[response]]==levels(data[[response]])[1])}

  if (is.null(subset)){subset <- rep(TRUE, nrow(data))}

  # Grow tree
  nodes <- growTemp(id=1L, depth=1L, data=data, response=response, treatment=treatment, Propensity=Propensity, subset=subset, search=search, method=method, split=split,
                    mtry=mtry, nsplit=nsplit, nsplit.random=nsplit.random, minsplit=minsplit, minbucket=minbucket,
                    maxdepth=maxdepth, a=a, scale.y=scale.y, trtlevels=trtlevels,response.type = response.type)

  # Compute terminal node number for each observation
  fitted <- fitted_node(nodes, data=data)
  thing <- ncol(data[Propensity])
  if(thing !=1) {
    daprop<-cbind(data[[treatment]],data[[Propensity]])
    ps<-apply(daprop,1,function(v){x<-v[1];return(v[x+1])})
  }
  ps<- data[[Propensity]]
  # Return rich constparty object
  ret <- party(nodes, data = data,
               fitted = data.frame("(fitted)" = fitted,
                                   "(response)" = data[[response]],
                                   "(treatment)" = data[[treatment]],
                                   "(propensity)"= ps,#data[[Propensity]],
                                   "(iptw)"= data[[Iptw]],
                                   check.names = FALSE),
               terms = terms(formula))
  as.constparty(ret)
  #as.simpleparty(ret)
}
