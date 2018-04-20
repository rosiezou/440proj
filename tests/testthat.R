library(testthat)
library(stat440pkg)

test_that("load data", {
  multiis.data <- multiis
  expect_that(multiis.data, is_a("data.frame"))
})

test_that("generate imputed responses",{
  imputed <- gen.imp.resp(multiis, num.iter = 1)
  expect_that(imputed, is_a("data.frame"))
})

test_that("generate latent variables",{
  grp.indicator <- sapply(names(multiis), FUN =
                            function(x){strsplit(x, split = "_")[[1]][2]})
  latent.vars <- gen.latent.vars(multiis, grp.indicator = grp.indicator, num.iter = 1)
  expect_that(latent.vars, is_a("data.frame"))
  expect_that(dim(latent.vars)[1], equals(333))
  expect_that(dim(latent.vars)[2], equals(3))
})


test_that("generate sets of latent variables",{
  grp.indicator <- sapply(names(multiis), FUN =
                            function(x){strsplit(x, split = "_")[[1]][2]})

  latent.datasets <- gen.latent.datasets(5, multiis, grp.indicator = grp.indicator, num.iter = 1)
  expect_that(latent.datasets, is_a("list"))
  expect_that(length(latent.datasets), equals(5))
  expect_that(dim(latent.datasets[[1]])[1], equals(333))
  expect_that(dim(latent.datasets[[1]])[2], equals(3))
})

test_that("pool analyses test",{
  grp.indicator <- sapply(names(multiis), FUN =
                            function(x){strsplit(x, split = "_")[[1]][2]})

  latent.datasets <- gen.latent.datasets(5, multiis, grp.indicator = grp.indicator, num.iter = 1)

  lm.pool <- pool.analyses(latent.datasets, cat~comp+int, lm)
  glm.pool <- pool.analyses(latent.datasets, cat~comp+int, glm)
  lm.pool.with.1.indepvar <- pool.analyses(latent.datasets, cat~comp, lm)
  glm.pool.with.1.indepvar <- pool.analyses(latent.datasets, cat~comp, glm)
  expect_that(lm.pool, is_a("list"))
  expect_that(glm.pool, is_a("list"))
  expect_that(lm.pool.with.1.indepvar, is_a("list"))
  expect_that(glm.pool.with.1.indepvar, is_a("list"))
  expect_that(length(lm.pool), equals(10))
  expect_that(length(glm.pool), equals(10))
  expect_that(length(lm.pool.with.1.indepvar), equals(10))
  expect_that(length(glm.pool.with.1.indepvar), equals(10))
  expect_that(dim(matrix(lm.pool$point.estimate))[1], equals(3))
  expect_that(dim(matrix(lm.pool.with.1.indepvar$point.estimate))[1], equals(2))
})
