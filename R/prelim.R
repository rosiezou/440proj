library(truncnorm)
library(MASS)

#' Obtain imputed continuous responses from ordinal responses
#'
#' Given a dataset of ordinal values (e.g. survey responses), obtain continuous imputed responses based
#' on quantiles of N(0,1), and using all other columns (similar to the MICE package).
#'
#' @importFrom MASS polr
#' @importFrom truncnorm rtruncnorm
#' @importFrom stats quantile qnorm simulate
#' @param data a data.frame containing ordinal variables (e.g. survey responses).
#' @param num.iter a number specifying the number of times to iterate the imputation (defaults to 20)
#' @return A data.frame containing the imputed responses.
#' @export
#' @examples
#' imputed <- gen.imp.resp(multiis, num.iter = 1)
#'
#' head(imputed)
gen.imp.resp <- function(data, num.iter = 5)
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
#' @return A data.frame containing the underlying latent variable estimates.
#' @export
#' @examples
#' # Create indicators (a label indicating which latent variable the question corresponds to)
#' grp.indicator <- sapply(names(multiis), FUN =
#'                          function(x){strsplit(x, split = "_")[[1]][2]})
#'
#' latent.vars <- gen.latent.vars(multiis, grp.indicator = grp.indicator, num.iter = 5)
#'
#' head(latent.vars)
gen.latent.vars <- function(data, grp.indicator, scores = "Bartlett", num.iter = 5)
{
  imp.resp.data <- gen.imp.resp(data, num.iter = num.iter)
  latent.labels <- unique(grp.indicator) # Get unique latent variable labels

  # Create as matrix first, then coerce to data.frame. This is because it's easiest
  # to make a matrix with our desired dimensions, and then treat it as a data.frame
  # later on.
  latent.vars.mtx <- matrix(nrow = nrow(imp.resp.data), ncol = length(latent.labels))
  latent.vars.df <- as.data.frame(latent.vars.mtx)
  names(latent.vars.df) <- latent.labels
  faresults <- vector(mode="list", length = length(latent.labels)) # Where to store results from factanal
  for(label in latent.labels)
  {
    latent.cols.indic <- grp.indicator == label # TRUE if column corresponds to label
    latent.cols.df <- imp.resp.data[latent.cols.indic] # sub-data.frame of only those columns corresp. to label
    fa <- factanal(latent.cols.df, factors = 1, scores = scores) # perform factor analysis
    latent.vars.df[label] <- fa$scores[,1] # Get scores *as vector* so we don't update column name
  }
  return(latent.vars.df)
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
#' # Create indicators (a label indicating which latent variable the question corresponds to)
#' grp.indicator <- sapply(names(multiis), FUN =
#'                          function(x){strsplit(x, split = "_")[[1]][2]})
#'
#' latent.datasets <- gen.latent.datasets(5, multiis, grp.indicator = grp.indicator, num.iter = 1)
#'
#' head(latent.datasets[[1]])
#'
#' head(latent.datasets[[2]])
gen.latent.datasets <- function(M, data, grp.indicator, scores = "Bartlett", num.iter = 5)
{
  latent.labels <- unique(grp.indicator) # labels or "names" of the different latent variables. e.g. c("cat", "int", "comp")
  num.vars <- length(latent.labels) # number of latent vars

  empty <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = num.vars)) # empty data.frame. Just initializing size here.
  names(empty) <- latent.labels # Name columns. Specifically each column will correspond to a latent variable
  datasets <- rep(list(empty), M) # Create M copies of the empty dataset in a list

  for(i in 1:M)
  {
    latent.vars <- gen.latent.vars(data = data, grp.indicator = grp.indicator,
                                   scores = scores, num.iter = num.iter) # Generate latent variables

    datasets[[i]] <- latent.vars # Store each generated latent dataset in a list
  }

  return(datasets)
}

#' Pool analyses results given M latent variable data sets, and estimate parameters
#'
#' @importFrom stats as.formula coef vcov pnorm lm
#' @param latent.datasets a (non-empty) list of lists returned by the gen.latent.vars function.
#' @param formula a valid formula in the form of response.var ~ independent.vars
#' @param method a regression function (e.g. lm). Currently only supports lm
#' @return A list of 10 elements consisting of: point estimate for parameter Q, within-imputaiton variance,
#' estimates of parameter Q obtained from M multiple imputations, estimates of variance, difference
#' between estimates of parameter Q and the final point estimate, between-imputation variance,
#' total variance, relative increase in variance due to nonresponse, fraction of missing information,
#' and results from hypothesis testing with H_0: beta_i = 0 for i = 1 ... number_of_parameters.
#' The mean and variance are calculated according to Rubin's Rules for Multiple Imputation. For more mathematical
#' details, please refer to page 5 of this UCLA paper
#' https://stats.idre.ucla.edu/wp-content/uploads/2016/02/multipleimputation.pdf \(section title
#' "Combining Inferences from Imputed Data Sets"\)
#' @export
#' @examples
#' # Create indicators (a label indicating which latent variable the question corresponds to)
#' grp.indicator <- sapply(names(multiis), FUN =
#'                          function(x){strsplit(x, split = "_")[[1]][2]})
#'
#' latent.datasets <- gen.latent.datasets(5, multiis, grp.indicator = grp.indicator, num.iter = 1)
#'
#' lm.pool <- pool.analyses(latent.datasets, cat~comp+int, lm)
#'
pool.analyses <- function(latent.datasets, formula, method = lm){
  ## m is the number of imputed data sets, requires at least 1 for pool analysis
  m <- length(latent.datasets)
  if (m < 1){
    stop("At least 1 analysis is needed for pooling")
  }

  ## A valid formula is required for pooled analysis. Currently only linear regression models
  ##    are supporte by this package. Formula must follow the format "dependent_var ~ independent_vars"
  if (is.null(formula)){
    stop("A valid formula is required.")
  }

  ## We designed the function so that in the future it could be easily extended to accept more methods.
  ##   So far only linear regression is supported.
  ## method cannot be null and must be lm
  if (is.null(method)){
    stop("A valid method is required.")
  }

  ## Apply the lm function to all M imputed data sets. Note that all data sets are fitted using
  ##   the same formula and lm function call
  fitted.objects <- lapply(latent.datasets,
                           FUN = function(x){
                             method(as.formula(formula), data = x)
                           }
  )

  ## Get the number and names of parameters
  k = length(fitted.objects[[1]]$coefficients)
  names <- names(coef(fitted.objects[[1]]))
  qhat <- matrix(NA, nrow = m, ncol = k, dimnames = list(seq_len(m),
                                                         names))
  u <- array(NA, dim = c(m, k, k), dimnames = list(seq_len(m),
                                                   names, names))

  ## This part of the code emulates the behaviour of pool() function in MICE package
  for (i in 1:m){
    fit <- fitted.objects[[i]] ## get fitted object
    qhat[i, ] <- coef(fit)     ## parameters for a single fitted object out of m
    ui <- vcov(fit)            ## variance-covariance matrix of the fitted object
    if (ncol(ui) != ncol(qhat))  ## double-check to make sure that we don't have dimensionality issues
      stop("Different number of parameters: coef(fit): ",
           ncol(qhat), ", vcov(fit): ", ncol(ui))
    u[i, , ] <- array(ui, dim = c(1, dim(ui)))
  }

  ## calcualte sample mean (i.e. point estimate) and sample variance using Rubin's rules
  qbar <- apply(qhat, 2, mean)          ## sample mean of parameters
  ubar <- apply(u, c(2, 3), mean)       ## within-imputation variance
  e <- qhat - matrix(qbar, nrow = m, ncol = k, byrow = TRUE)   ## difference betw. each parameter estimate and mean
  b <- (t(e) %*% e)/(m - 1)             ## between-imputation variance
  t <- ubar + (1 + 1/m) * b             ## total variance
  r <- (1 + 1/m) * diag(b/ubar)         ## relative increase in variance due to non-response
  lambda <- (1 + 1/m) * diag(b/t)       ## fraction of missing information

  results <- matrix(NA, nrow = k, ncol = 3) ## package the results to be displayed
  rownames(results) <- names
  colnames(results) <- c("Estimate", "Std. Error", "p value")
  p <- matrix(NA, nrow = k, ncol = 1)

  for (i in 1:k){
    se = sqrt(t[i,i])/sqrt(m)          ## calculate standard error
    pval <- 2*pnorm(abs(qbar[i]/se), lower.tail=FALSE)  ## calculate p values
    p[i,1] <- pval
  }

  results[, 1] = t(qbar)
  results[, 2] = t(sqrt(diag(t))/sqrt(m))
  results[, 3] = p

  results <- signif(results, 3)

  rval <- list(point.estimate = qbar, within.imputation.variance = ubar,
               multiple.imputation.estimates = qhat, variance.estimates = u,
               differences.from.point.estimate = e, between.imputation.variance = b,
               total.variance = t, relative.increase.in.variance.due.to.nonresponse = r,
               fraction.of.missing.info = lambda, hypothesis.test = results)

  return(rval)
}
