load_wood <- function() {
  read.csv(system.file("extdata/wood.csv", package=.packageName,
                       mustWork=TRUE),
           stringsAsFactors=FALSE)
}
