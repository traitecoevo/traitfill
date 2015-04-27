# traitfill

[![Build Status](https://travis-ci.org/richfitz/traitfill.png?branch=master)](https://travis-ci.org/richfitz/traitfill)
[![Coverage Status](https://coveralls.io/repos/richfitz/traitfill/badge.svg?branch=master)](https://coveralls.io/r/richfitz/traitfill?branch=master)

This R package is for leveraging taxonomic information to impute the values of a binary trait for species with missing data. The approach is described in [FitzJohn et al. 2014 (Journal of Ecology)](http://onlinelibrary.wiley.com/doi/10.1111/1365-2745.12260/full) and originally implemented [here](https://github.com/richfitz/wood).

The easiest way to install the package is with `devtools`
``` 
## install.packages("devtools")
devtools::install_github("richfitz/traitfill")
```

The primary function of the package is also called `traitfill`. To demonstrate how it works, we will use the data from our woodiness study
```
library(traitfill)
wood <- load_wood()
res <- traitfill(wood, nrep=50, with_replacement=FALSE, names=c("H", "W"))
res
```

As described in FitzJohn et al. 2014, we consider two different sampling distributions for the traits:

1. Sampling with replacement (`with_replacement = TRUE` in the `traitfill` function) in which trait values for missing taxa are drawn from a **binomial distribution**. If all the records in a genus are of one type, then this implies that the rest of the unsampled species also share this trait.

2. Sampling without replacement (`with_replacement = FALSE` in the `traitfill` function) in which trait values for missing taxa are drawn from a **hypergeometric distribution**. Even if all trait records for a genus are of one type, we consider the possibility that there may be a species of the other type that we simply have not sampled. The probablility of this descreases as the proportion of taxa that are sampled increases.

