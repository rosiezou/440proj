library(truncnorm)
library(MASS)

#' Obtain imputed continuous responses from ordinal responses
#'
#' Given a dataset of ordinal values (e.g. survey responses), obtain continuous imputed responses based
#' on quantiles of N(0,1), and using all other columns (similar to the MICE package).
#'
#' @importFrom MASS polr
#' @importFrom truncnorm rtruncnorm
#' @importFrom stats quantile
#' @param data a data.frame containing ordinal variables (e.g. survey responses).
#' @param num.iter a number specifying the number of times to iterate the imputation (defaults to 20)
#' @return A data.frame containing the imputed responses.
#' @export
#' @examples
#' multiis <- read.csv("data/MULTIIS.csv")
#' imputed <- gen.imp.resp(multiis)
#'
#' head(imputed)
gen.imp.resp <- function(data, num.iter = 20)
{
  data.interp <- data # For initial interpolation

  for(q in names(data))
  {
    responses <- data[[q]]
    min.resp <- min(responses)
    max.resp <- max(responses)

    probs <- cumsum(colSums(sapply(min.resp:max.resp,
                                   FUN = function(i) {responses == i})) / nrow(data))

    quantiles <- c(-Inf, qnorm(probs)) # Get all quantiles at the probs, including -Inf

    lower <- quantiles[responses] # Get lower bound of bins for responses
    upper <- quantiles[responses + 1] # Same but get upper bound

    # Generate from truncated normal
    sim.truncnorm <- rtruncnorm(n = 1, a = lower, b = upper, mean = 0, sd = 1)

    data.interp[q] <- sim.truncnorm
  }

  data.analysis <- data.interp

  for (i in 1:num.iter) # number of updates
  {
    for(q in names(data.analysis)) # For each column
    {
      data.analysis[q] <- data[[q]]

      # Construct regression formula with q as the response
      regression.formula <- as.formula(paste("as.ordered(", q, ")", "~ ."))
      inhs.polr <- polr(regression.formula, data = data.analysis, method = "probit")
      quantiles <- c(-Inf, inhs.polr$zeta, Inf) # Get quantiles that include -Inf and Inf

      sim.bins <- as.numeric(simulate(inhs.polr)[[1]]) # Simulate new bins

      lower <- quantiles[sim.bins] # Get lower bounds of bins in sim.bins
      upper <- quantiles[sim.bins + 1] # Get upper bounds of bins in sim.bins

      # Generate from truncated normal
      sim.truncnorm <- rtruncnorm(n = 1, a = lower, b = upper, mean = 0, sd = 1)

      data.analysis[q] <- sim.truncnorm # Update column
    }
  }

  return(data.analysis)
}

#' Generate latent underlying variables given imputed question responses
#'
#' Given ordinal responses to survey questions, we generate the
#' underlying latent variables that correspond to those questions
#'
#' @importFrom stats factanal
#' @param data a (non-empty) numeric vector of data values.
#' @param grp.indicator a (non-empty) numeric vector of data values.
#' @param scores type of score to use.
#' @param num.iter number of iterations to use in imputation step
#' @return A list l of length 2, where l[[1]] is a data.frame containing the underlying
#' latent variable estimates, and l[[2]] is a vector of j factor analyses objects
#' (j is the total number of latent variables).
#' @export
#' @examples
#' setwd("~/GitProjects/440proj") # Set to project folder
#' multiis <- read.csv("data/MULTIIS.csv")
#' imputed <- gen.imp.resp(multiis)
#' # Create indicators (a label indicating which latent variable the question corresponds to)
#' grp.indicator <- sapply(names(imputed), FUN =
#'                          function(x){strsplit(x, split = "_")[[1]][2]})
#'
#' latent.vars <- gen.latent.vars(imputed, grp.indicator = grp.indicator)
#'
#' head(latent.vars)
gen.latent.vars <- function(data, grp.indicator, scores = "Bartlett", num.iter = 20)
{
  imp.resp.data <- gen.imp.resp(data, num.iter = num.iter)
  latent.labels <- unique(grp.indicator) # Get unique latent variable labels

  # Create as matrix first, then coerce to data.frame. This is because it's easiest
  # to make a matrix with our desired dimensions, and then treat it as a data.frame
  # later on.
  latent.vars.mtx <- matrix(nrow = nrow(imp.resp.data), ncol = length(latent.labels))
  latent.vars.df <- as.data.frame(latent.vars.mtx)
  names(latent.vars.df) <- latent.labels
  faresults <- vector(mode="list", length = length(latent.labels))
  index = 1
  for(label in latent.labels)
  {
    latent.cols.indic <- grp.indicator == label # TRUE if column corresponds to label
    latent.cols.df <- imp.resp.data[latent.cols.indic] # sub-data.frame of only those columns corresp. to label

    fa <- factanal(latent.cols.df, factors = 1, scores = scores) # perform factor analysis

    faresults[[index]] <- fa
    index <- index + 1
    latent.vars.df[label] <- fa$scores[,1] # Get scores *as vector* so we don't update column name
  }

  return(list(latent.vars.df, faresults))
}

#' Generate multiple datasets of latent underlying variables
#'
#' Given ordinal responses to survey questions, we generate multiple
#' underlying latent variable datasets
#'
#' @param M number of datasets to construct.
#' @param data a (non-empty) numeric vector of data values.
#' @param grp.indicator a (non-empty) numeric vector of data values.
#' @param scores type of score to use.
#' @param num.iter number of iterations to use in imputation step
#' @return A list containing datasets with latent underlying variables
#' @export
#' @examples
#' setwd("~/GitProjects/440proj") # Set to project folder
#' multiis <- read.csv("data/MULTIIS.csv")
#'
#' # Create indicators (a label indicating which latent variable the question corresponds to)
#' grp.indicator <- sapply(names(multiis), FUN =
#'                          function(x){strsplit(x, split = "_")[[1]][2]})
#'
#' latent.datasets <- gen.latent.datasets(5, multiis, grp.indicator = grp.indicator, num.iter = 5)
#'
#' head(latent.datasets[[1]])
#'
#' head(latent.datasets[[2]])
gen.latent.datasets <- function(M, data, grp.indicator, scores = "Bartlett", num.iter = 20)
{
  latent.labels <- unique(grp.indicator)
  num.vars <- length(latent.labels)

  empty <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = num.vars))
  names(empty) <- latent.labels
  datasets <- rep(list(empty), M)

  for(i in 1:M)
  {
    latent.vars <- gen.latent.vars(data = data, grp.indicator = grp.indicator,
                                   scores = scores, num.iter = num.iter)[[1]]

    datasets[[i]] <- latent.vars
  }

  return(datasets)
}

#' Pool analyses results given M latent variable data sets, and estimate parameters
#'
#' @param latent.datasets a (non-empty) list of lists returned by the gen.latent.vars function.
#' @return A list of lists of 5 elements consisting of: within-imputation variance, between-imputation variance, total variance, number of analyses, and a vector of final estimated parameters
#' @export
#' @examples
#'
pool.analyses <- function(latent.datasets){
  M <- length(latent.datasets)
  if (M < 1){
    stop("At least 1 analysis is needed for pooling")
  }

  n_latent_vars <- length(latent.datasets[[1]][[2]]) ## get the total number of latent variables

  pool.out <- vector(mode="list", length = n_latent_vars) ## final result will be a list of length n_latent_vars

  for (i in 1:n_latent_vars){

    final_loading_matrix <- c() ## used for calculating between-imputation variance for latent variable i

    sample_variances <- 0 ## used for calculating within-imputation variance for latent variable i

    for (j in 1:M){
      imputed.data <- latent.datasets[[j]][[1]]
      fa.obj <- latent.datasets[[j]][[2]][[i]]
      row <- t(fa.obj$loadings)
      final_loading_matrix <- rbind(final_loading_matrix, row)
      sample_variances <- sample_variances + var(imputed.data[,i])
    }

    within_imputation_variance <- sample_variances / M

    between_imputation_variance <- 0

    mean_parameters <- final_loading_matrix[1,]

    for (k in 1:M){
      if (k > 1){
        mean_parameters <- mean_parameters + final_loading_matrix[k,]
      }
    }

    mean_parameters <- mean_parameters / M

    for (l in 1:M){
      vec <- final_loading_matrix[l,] - mean_parameters
      between_imputation_variance <- between_imputation_variance + (t(vec) %*% vec)[1,1]
    }

    between_imputation_variance <- between_imputation_variance / (M - 1)

    total_variance <- within_imputation_variance + between_imputation_variance +
      between_imputation_variance/M ## total variance as calculated per Rubin's rule

    ## now pack all key values into one list and then add it to pool.out
    pool.out[[i]] <- list(within_imputation_variance, between_imputation_variance,
                          total_variance, M, mean_parameters)
  }
  return(pool.out)
}

#' Perform a correlation test on multiply imputed datasets
#'
#' Given imputed datasets, perform a correlation test between specified columns
#'
#' @param datasets a list of data.frames
#' @param indices a vector of length two specifying the index of the columns on which to test correlation
#' @param grp.indicator a vector indicating which underlying group the columns in each dataset correpsonds to
#' @param alternative the alternative hypothesis
#' @param method the method to compute correlation. Currently only "pearson" is supported.
#' @return A list containing the mean correlation and a p-value from the hypothesis test
#' @export
#' @examples
#' setwd("~/GitProjects/440proj") # Set to project folder
#' multiis <- read.csv("data/MULTIIS.csv")
#'
#' # Create indicators (a label indicating which latent variable the question corresponds to)
#' grp.indicator <- sapply(names(multiis), FUN =
#'                          function(x){strsplit(x, split = "_")[[1]][2]})
#'
#' latent.datasets <- gen.latent.datasets(5, multiis, grp.indicator = grp.indicator, num.iter = 5)
#'
#' pooled.cor.test(latent.datasets, indices = c(1,3))
pooled.cor.test <- function(datasets, indices = c(1, 2), alternative = "two.sided", method = "pearson")
{
  num.sets <- length(datasets) # number of datasets
  num.row <- nrow(datasets[[1]]) # number of rows in each dataset
  corr <- numeric(num.sets)
  piv.val <- numeric(num.sets) # pivotal value

  for(i in 1:num.sets)
  {
    dataset <- datasets[[i]]
    if(method == "pearson")
    {
      r.value <- cor(dataset[indices[1]], dataset[indices[2]]) # sample correlation
      corr[i] <- r.value
      piv.val[i] <- r.value * sqrt(num.row - 2) / sqrt(1 - r.value^2)
    }
  }

  piv.var.between <- sum((piv.val - mean(piv.val))^2) / (num.sets - 1)
  piv.var.total <- 1 + piv.var.between + piv.var.between / num.sets # var.within is 1 since distribution is t, so asymptotically N(0,1)

  p.val <- 2 * pnorm(abs(mean(piv.val) / sqrt(piv.var.total)), lower.tail = FALSE)

  return(list(mean.corr = mean(corr),
              mean.piv.val = mean(piv.val),
              piv.var.total = piv.var.total,
              piv.var.between = piv.var.between,
              p.val = p.val))
}
