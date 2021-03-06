vlapply <- function(X, FUN, ...) {
  vapply(X, FUN, logical(1), ...)
}
viapply <- function(X, FUN, ...) {
  vapply(X, FUN, integer(1), ...)
}
vnapply <- function(X, FUN, ...) {
  vapply(X, FUN, numeric(1), ...)
}
vcapply <- function(X, FUN, ...) {
  vapply(X, FUN, character(1), ...)
}

rename <- function(x, from, to) {
  i <- match(from, names(x))
  if (any(is.na(i))) {
    stop("Did not find names ", paste(from[is.na(i)], collapse=", "))
  }
  names(x)[i] <- to
  x
}

pastec <- function(...) {
  paste(..., collapse=", ")
}
