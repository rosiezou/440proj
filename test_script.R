library(stat440pkg)

grp.indicator <- sapply(names(multiis), FUN =
                          function(x){strsplit(x, split = "_")[[1]][2]})

latent.datasets <- gen.latent.datasets(100, multiis, grp.indicator = grp.indicator, num.iter = 5)
pooled.add1 <- pool.analyses(latent.datasets, cat~comp + int, lm)
pooled.add2 <- pool.analyses(latent.datasets, comp~cat + int, lm)
pooled.add3 <- pool.analyses(latent.datasets, int~comp + cat, lm)


library(Rcmdr)
add <- function(x) Reduce("+", x)
averaged <- add(latent.datasets)/100
scatter3d(averaged$comp, averaged$cat, averaged$int)