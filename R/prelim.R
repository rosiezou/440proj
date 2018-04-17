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
#' @return A list l of length 2, where l[[1]] is a data.frame containing the underlying
#' latent variable estimates, and l[[2]] is a vector of j factor analyses objects
#' (j is the total number of latent variables).
#' @export
#' @examples
#' # Create indicators (a label indicating which latent variable the question corresponds to)
#' grp.indicator <- sapply(names(multiis), FUN =
#'                          function(x){strsplit(x, split = "_")[[1]][2]})
#'
#' latent.vars <- gen.latent.vars(multiis, grp.indicator = grp.indicator)
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
  latent.labels <- unique(grp.indicator)
  num.vars <- length(latent.labels)

  empty <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = num.vars))
  names(empty) <- latent.labels
  datasets <- rep(list(empty), M)

  for(i in 1:M)
  {
    latent.vars <- gen.latent.vars(data = data, grp.indicator = grp.indicator,
                                   scores = scores, num.iter = num.iter)

    datasets[[i]] <- latent.vars
  }

  return(datasets)
}

#' Pool analyses results given M latent variable data sets, and estimate parameters
#'
#' @importFrom stats as.formula coef vcov
#' @param latent.datasets a (non-empty) list of lists returned by the gen.latent.vars function.
#' @param formula a valid formula in the form of response.var ~ independent.vars
#' @param method a regression function (e.g. lm). Currently only supports lm and glm
#' @return A list 9 elements consisting of: point estimate for parameter Q, within-imputaiton variance,
#' estimates of parameter Q obtained from M multiple imputations, estimates of variance, difference
#' between estimates of parameter Q and the final point estimate, between-imputation variance,
#' total variance, relative increase in variance due to nonresponse, and fraction of missing information.
#' These values are calculated according to Rubin's Rules for Multiple Imputation. For more mathematical
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
#' glm.pool <- pool.analyses(latent.datasets, cat~comp+int, glm)
#'
pool.analyses <- function(latent.datasets, formula, method){
  m <- length(latent.datasets)
  if (m < 1){
    stop("At least 1 analysis is needed for pooling")
  }

  if (is.null(formula)){
    stop("A valid formula is required.")
  }

  if (is.null(method)){
    stop("A valid method is required.")
  }

  ## currently only supports lm and glm(family = "gaussian")
  fitted.objects <- lapply(latent.datasets,
                           FUN = function(x){
                             method(as.formula(formula), data = x)
                           }
  )

  k = length(latent.datasets[[1]][1,])
  names <- names(coef(fitted.objects[[1]]))
  qhat <- matrix(NA, nrow = m, ncol = k, dimnames = list(seq_len(m),
                                                         names))
  u <- array(NA, dim = c(m, k, k), dimnames = list(seq_len(m),
                                                   names, names))

  for (i in 1:m){
    fit <- fitted.objects[[i]]
    qhat[i, ] <- coef(fit)
    ui <- vcov(fit)
    if (ncol(ui) != ncol(qhat))
      stop("Different number of parameters: coef(fit): ",
           ncol(qhat), ", vcov(fit): ", ncol(ui))
    u[i, , ] <- array(ui, dim = c(1, dim(ui)))
  }
  qbar <- apply(qhat, 2, mean)
  ubar <- apply(u, c(2, 3), mean)
  e <- qhat - matrix(qbar, nrow = m, ncol = k, byrow = TRUE)
  b <- (t(e) %*% e)/(m - 1)
  t <- ubar + (1 + 1/m) * b
  r <- (1 + 1/m) * diag(b/ubar)
  lambda <- (1 + 1/m) * diag(b/t)
  list(point.estimate = qbar, within.imputation.variance = ubar,
       multiple.imputation.estimates = qhat, variance.estimates = u,
       differences.from.point.estimate = e, between.imputation.variance = b,
       total.variance = t, relative.increase.in.variance.due.to.nonresponse = r,
       fraction.of.missing.info = lambda)
}

#' Perform a correlation test on multiply imputed datasets
#'
#' Given imputed datasets, perform a correlation test between specified columns
#'
#' @importFrom stats pnorm cor
#' @param datasets a list of data.frames
#' @param indices a vector of length two containing the indices of columns to test (either column name or number)
#' @param alternative the alternative hypothesis
#' @param method the method to compute correlation. Either "pearson" or "spearman"
#' @return A list containing the mean correlation and a p-value from the hypothesis test
#' @export
#' @examples
#' # Create indicators (a label indicating which latent variable the question corresponds to)
#' grp.indicator <- sapply(names(multiis), FUN =
#'                          function(x){strsplit(x, split = "_")[[1]][2]})
#'
#' latent.datasets <- gen.latent.datasets(5, multiis, grp.indicator = grp.indicator, num.iter = 1)
#'
#' pooled.cor.test(latent.datasets, indices = c("cat","int"), method = "pearson")
#' pooled.cor.test(latent.datasets, indices = c("cat", "int"), method = "spearman")
pooled.cor.test <- function(datasets, indices = c(1, 2), alternative = "two.sided", method = "pearson")
{
  num.sets <- length(datasets) # number of datasets
  num.row <- nrow(datasets[[1]]) # number of rows in each dataset
  corr <- numeric(num.sets)
  piv.val <- numeric(num.sets) # pivotal value

  for(i in 1:num.sets)
  {
    dataset <- datasets[[i]]

    r.value <- cor(dataset[indices[1]], dataset[indices[2]], method = method) # sample correlation
    corr[i] <- r.value
    piv.val[i] <- r.value * sqrt(num.row - 2) / sqrt(1 - r.value^2)
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
