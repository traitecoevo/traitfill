##' Fill missing binary trait data values based on taxonomic structure.
##'
##' This is a \emph{direct} port from richfitz/wood with nothing other
##' than stylistic changes.
##'
##' It assumes that the taxonomy is structured into order / family /
##' genus.  Trait data are assumed in columns called "n0" and "n1",
##' unless alternative names are passed through in \code{names}.  A
##' column called "N" must be present that has the estimate of the
##' number of species within the genus.
##'
##' So \code{dat} must have columns:
##' \itemize{
##'   \item{order}
##'   \item{family}
##'   \item{genus}
##'   \item{n1 (unless renamed by \code{names})}
##'   \item{n0 (unless renamed by \code{names})}
##'   \item{N (must be at least n0 + n1)}
##' }
##'
##' The function returns a big list with elements:
##' \itemize{
##'   \item{genus: estimates at the level of genus}
##'   \item{family: estimates at the level of family}
##'   \item{order: estimates at the level of order}
##'   \item{overall: overall estimates}
##' }
##'
##' each element is a data.frame with columns
##' \itemize{
##'   \item{n0, n1: observed numbers of species in each state,
##'     possibly renamed}
##'   \item{N: number of species at that taxonomic level}
##'   \item{K: number of species with known states}
##'   \item{mean, lower, upper: mean and lower/upper bounds on the
##'     estimated number of species in state 1}
##'   \item{p_mean, p_lower, p_upper: mean and lower/upper bounds on the
##'     estimated fraction of species in state 1}
##' }
##' 
##' @title Fill missing trait values
##' @param dat A data set (see Details)
##' @param nrep Number of simulation runs to do.
##' @param with_replacement Do the sampling with or without
##' replacement (default is TRUE)
##' @param names Names for columns in \code{dat} corresponding to
##' state 0 and state1.
##' @export
##' @examples
##' ## Data from "How much of the world is woody?"
##' wood <- load_wood()
##' ## Fill in missing woodiness data; this says that the "H" column
##' ## is state 0 and the "W" column is state 1, so that the generated
##' ## percentages are "percentage woody" and not "percentage herbacious.
##' res <- traitfill(wood, 50, names=c("H", "W"))
traitfill <- function(dat, nrep, with_replacement=TRUE,
                      names=c("n0", "n1")) {
  f <- function(level) {
    tmp <- lapply(res, "[[", level)
    ret <- do.call("rbind", tmp, quote=TRUE)
    rownames(ret) <- NULL
    ret[c("p_mean", "p_lower", "p_upper")] <-
      ret[c("mean", "lower", "upper")] / ret[["N"]]
    ret <- rename(ret, traitfill_names(), names)
    ret
  }

  dat <- traitfill_prepare(dat, names)

  res <- lapply(split(dat, dat$order), sim, nrep, with_replacement)

  total <- rowSums(sapply(res, attr, "total"))
  overall <- summarise(total)
  overall_p <- overall / sum(dat$N)

  overall <- as.data.frame(c(n0=sum(dat$n0),
                             n1=sum(dat$n1),
                             N=sum(dat$N),
                             K=sum(dat$n0, dat$n1),
                             as.list(summarise(total))))
  overall[c("p_mean", "p_lower", "p_upper")] <-
    overall[c("mean", "lower", "upper")] / sum(overall$N)
  overall <- rename(overall, traitfill_names(), names)

  ret <- list(genus=f("genus"),
       family=f("family"),
       order=f("order"),
       overall=overall)
}

sim <- function(x, nrep, with_replacement, names) {
  ## First, focus on cases where we have a valid estimate of the
  ## fraction of species that are woody (i.e., at least one known
  ## species).
  ok <- x$K > 0

  w <- matrix(NA, nrow(x), nrep)

  ## A: genera with any known species
  if (with_replacement) {
    w[ok, ] <- x$n1[ok] + rbinom(sum(ok), x$N[ok] - x$K[ok], x$n1[ok] / x$K[ok])
  } else {
    w[ok, ] <- t(sapply(which(ok),
                        function(i) rhyper2(nrep, x$n0[i], x$n1[i], x$N[i])))
  }

  ## B: genera with no known species
  n_unknown <- sum(!ok)
  g <- function(y) {
    rbinom(n_unknown, x$N[!ok], quantile(y, runif(n_unknown)))
  }
  w[!ok, ] <- apply(w[ok, , drop=FALSE] / x$N[ok], 2, g)
  
  rownames(w) <- x$genus

  summarise_sim(w, x[c("order", "family", "genus",
                       "n1", "n0", "N", "K")])
}

## This collects up the results at different taxonomic levels.
summarise <- function(x, p=1/20) {
  structure(c(mean(x), quantile(x, c(p / 2, 1 - p/2))),
            names=c("mean", "lower", "upper"))
}

summarise_sim <- function(w, info) {
  order <- info$order[[1]]

  info_cols <- c("n1", "n0", "N", "K")

  ## Genus is easy;
  w_g <- cbind(info, t(apply(w, 1, summarise)))

  ## Family is a pain:
  w_f <- do.call(rbind,
                 lapply(split(as.data.frame(w), info$family), colSums))
  w_f <- t(apply(w_f, 1, summarise))
  w_f <- data.frame(order=order,
                    aggregate(info[info_cols], info["family"], sum),
                    w_f, stringsAsFactors=TRUE) # really?
  rownames(w_f) <- NULL
  
  ## Order is easy; we are guaranteed to have just one order here, so:
  w_o <- data.frame(order=order,
                    as.data.frame(t(colSums(info[info_cols]))),
                    t(summarise(colSums(w))), stringsAsFactors=FALSE)

  ret <- list(genus=w_g, family=w_f, order=w_o)
  attr(ret, "total") <- colSums(w)
  ret
}

rhyper2 <- function(nn, s0, s1, xn, fraction=FALSE) {
  x1 <- seq(s1, xn - s0)
  x0 <- xn - x1
  p1 <- dhyper(s1, x1, x0, s0+s1)
  p1 <- p1 / sum(p1)
  x1[sample(length(p1), nn, TRUE, p1)]
}

traitfill_prepare <- function(dat, names) {
  taxon <- c("order", "family", "genus")
  msg_taxon <- setdiff(c("order", "family", "genus"), names(dat))
  if (length(msg_taxon) > 0L) {
    stop("Missing taxonomic columns: ", pastec(msg_taxon))
  }
  dat[taxon] <- lapply(dat[taxon], as.character)
  # TODO: this *should* be correct but there is an issue with the
  # woodiness data needs fixing first.
  # msg <- any(is.na(dat[taxon])) || any(dat[taxon] == "")
  msg <- any(dat[taxon] == "", na.rm=TRUE)
  if (msg) {
    stop("Detected missing taxonomic information")
  }

  msg_data <- setdiff(c(names, "N"), names(dat))
  if (length(msg_data) > 0L) {
    stop("Missing data columns: ", pastec(msg_data))
  }

  ## rename the data columns to standard names to work with
  ## throughout:
  dat <- rename(dat, names, traitfill_names())

  dat$K <- dat$n0 + dat$n1

  i <- dat$N < dat$K
  if (any(i)) {
    stop("Some rows have less species than data: ", pastec(i))
  }

  dat
}

traitfill_names <- function() {
  c("n0", "n1")
}
