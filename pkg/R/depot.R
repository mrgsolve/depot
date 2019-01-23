##' Model depot location
##' @importFrom mrgsolve mread
##' @export
depot <- function(x = NULL, ...) {
  if(is.null(x)) return(system.file("models", package = "depot"))
  mrgsolve::mread(x, depot(),...)
}
