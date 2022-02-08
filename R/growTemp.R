### Grow tree by using partition() function several times in a recursive loop ###
growTemp <- function(id=1L, depth=1L, data, response, treatment, Propensity, subset, search, method, split,
                     mtry, nsplit, nsplit.random, minsplit, minbucket, maxdepth,
                     a, scale.y, trtlevels,response.type){
  if (depth > maxdepth) {return(partynode(id=id))}
  y <- data[[response]]
  trt<-data[[treatment]]
  propensity<-data[[Propensity]]
  varSelected <- sort(sample.int(ncol(data)-4, mtry))
  vars <- data[varSelected]
  colnames(vars) <- varSelected #Have columns represent varid

  sp <- partition(vars=vars, y=y,  subset=subset,trt=trt,propensity=propensity,
                  search=search, method=method, split=split, nsplit=nsplit, nsplit.random=nsplit.random,
                  minsplit=minsplit, minbucket=minbucket, a=a, scale.y=scale.y,
                  trtlevels=trtlevels,response.type = response.type)
  #useSearch=useSearch, useOptim=useOptim,
  if (is.null(sp)) {return(partynode(id=id))}

  # Split the data
  kidids <- kidids_split(sp, data=data)
  depth <- depth + 1
  #print(max(kidids, na.rm=TRUE))
  kids <- vector(mode="list", length=max(kidids, na.rm=TRUE))
  for (kidid in seq_along(kids)) {
    s <- subset # subset is the previous loops s
    s[kidids != kidid] <- FALSE
    # Node ID
    if (kidid > 1) {myid <- max(nodeids(kids[[kidid-1]]))
    } else {myid <- id}
    # Start recursion on this daugther node
    kids[[kidid]] <- growTemp(id=as.integer(myid+1), depth=depth, data=data, response=response, treatment=treatment, Propensity=Propensity,
                              subset=s, search=search, method=method, split=split, mtry=mtry, nsplit=nsplit, nsplit.random=nsplit.random,
                              minsplit=minsplit, minbucket=minbucket, maxdepth=maxdepth,
                              a=a, scale.y=scale.y, trtlevels=trtlevels,response.type = response.type)
  }
  #print(sapply(kids, class))
  #print(length(kids))
  #print(class(kids))
  return(partynode(id=as.integer(id), split=sp, kids=kids,
                   #info=list(stats=max(info_split(sp)$stats, na.rm=TRUE))))
                   info = depth))
}
