library(stat440pkg)
multiis <- read.csv("data/MULTIIS.csv")
M = 100
all_imputations <- vector(mode = "list", length = M)
for (i in 1:M){
  imputed <- gen.imp.resp(multiis)
  grp.indicator <- sapply(names(imputed), FUN =
                            function(x){strsplit(x, split = "_")[[1]][2]})
  latent.vars <- gen.latent.vars(imputed, grp.indicator = grp.indicator)
  all_imputations[[i]] <- latent.vars
}

pooled_results <- pool.analyses(M, all_imputations)