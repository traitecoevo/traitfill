context("traitfill")

test_that("woodiness", {
  wood <- load_wood()
  set.seed(1)
  res <- traitfill(wood, 2, names=c("H", "W"))

  expect_that(names(res), equals(c("genus", "family", "order",
                                   "overall")))
  expect_that("H" %in% names(res$genus), is_true())
  expect_that("W" %in% names(res$genus), is_true())

  r <- c(0.44, 0.47)
  
  expect_that(res$overall$p_mean, is_more_than(r[[1]]))
  expect_that(res$overall$p_mean, is_less_than(r[[2]]))

  set.seed(1)
  res2 <- traitfill(wood, 2, names=c("W", "H"))

  expect_that(res2$overall$p_mean, is_more_than(1 - r[[2]]))
  expect_that(res2$overall$p_mean, is_less_than(1 - r[[1]]))
})
