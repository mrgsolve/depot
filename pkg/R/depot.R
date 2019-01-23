##' Access the model repo
##'
##' Pass in the name of a model to load the model; or call with
##' no arguments to get the path to those models (like [mrgsolve::modlib]).
##'
##'
##' @importFrom mrgsolve mread
##'
##' @param x character name of a model in the depot
##' @param ... passed to [mrgsolve::mread]
##'
##' @examples
##'
##' \dontrun{
##'
##' mod <- depot("gcsf")
##'
##' mod <- mread("moxi", depot())
##' }
##'
##' @md
##' @export
depot <- function(x = NULL, ...) {
  if(is.null(x)) return(system.file("models", package = "depot"))
  mrgsolve::mread(x, depot(),...)
}
