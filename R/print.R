#' @export
print.CERFIT <- function(x,...){
  cat(paste("Numer of Trees:",length(x),"\n"))
  cat(paste("Treatment Type:",x[[1]]$trt.type))
}
# CapStr <- function(y) {
#   c <- strsplit(y, " ")[[1]]
#   paste(toupper(substring(c, 1,1)), substring(c, 2),
#         sep="", collapse=" ")
# }
