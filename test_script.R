library(stat440pkg)
M <- 50
grp.indicator <- sapply(names(multiis), FUN =
                          function(x){strsplit(x, split = "_")[[1]][2]})
latent.datasets <- gen.latent.datasets(M, multiis, grp.indicator = grp.indicator, num.iter = 5)
pooled.add1 <- pool.analyses(latent.datasets, cat~comp + int, lm)
pooled.add2 <- pool.analyses(latent.datasets, comp~cat + int, lm)
pooled.add3 <- pool.analyses(latent.datasets, int~comp + cat, lm)


library(scatterplot3d)
add <- function(x) Reduce("+", x)
averaged <- add(latent.datasets)/M
fit <- lm(int~comp + cat, data = averaged)
scplot <- scatterplot3d(averaged$comp, averaged$cat, averaged$int,
              main="3D Scatterplot of Latent Variables\n with Regression Plane for Int ~  Comp + Cat", angle = 120,
              xlab = "compartmentalization", ylab = "categorization", zlab = "integration",
              col.grid = "lightgrey", pch = 19, color = "lightblue")

scplot$plane3d(fit, lty = "dotted")
orig <- scplot$xyz.convert(averaged$comp, averaged$cat, averaged$int)
plane <- scplot$xyz.convert(averaged$comp, averaged$cat, fitted(fit))
i.negpos <- 1 + (resid(fit) > 0)
segments(orig$x, orig$y, plane$x, plane$y,
         col = c("blue", "red")[i.negpos], lty = (2:1)[i.negpos])

pairs(averaged)
plot(fit)

# ggplot2 pairs plot
library(ggplot2)
library(GGally)

ggpairs(averaged)
