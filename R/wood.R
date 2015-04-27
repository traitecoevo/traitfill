##' Data from the woodiness paper
##' @title Data from the woodiness paper
##' @export
##' @examples
##' head(load_wood())
load_wood <- function() {
  read.csv(system.file("extdata/wood.csv", package=.packageName,
                       mustWork=TRUE),
           stringsAsFactors=FALSE)
}
