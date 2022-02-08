#' @export
print.CERFIT <- function(x,...){
  cat(paste("Numer of Trees:",length(x$randFor),"\n"))
  cat(paste("Treatment Type:",x$trt.type, "\n"))
  cat(paste("Response Type:",x$response.type))
}
# CapStr <- function(y) {
#   c <- strsplit(y, " ")[[1]]
#   paste(toupper(substring(c, 1,1)), substring(c, 2),
#         sep="", collapse=" ")
# }
