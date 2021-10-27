print.CERFIT <- function(cerfit){
  cat(paste("Numer of Trees:",length(cerfit),"\n"))
  cat(paste("Treatment Type:",CapStr(cerfit[[1]]$trt.type)))
}
CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}
