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
#' @importFrom data.table setnames
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

#' Generate latent underlying variables given latent question responses
#'
#' Given imputed latent variable responses to survey questions, we generate the
#' underlying latent variables that correspond to those questions
#'
#' @importFrom stats factanal
#' @param data a (non-empty) numeric vector of data values.
#' @param grp.indicator a (non-empty) numeric vector of data values.
#' @param scores number of bootstrap resamples to perform.
#' @return A data.frame containing the underlying latent variable estimates.
#' @export
#' @examples
#' # Create indicators (a label indicating which latent variable the question corresponds to)
#' grp.indicator <- sapply(names(imputed), FUN =
#'                          function(x){strsplit(x, split = "_")[[1]][2]})
#'
#' latent.vars <- gen.latent.vars(imputed, grp.indicator = grp.indicator)
#'
#' head(latent.vars)
gen.latent.vars <- function(data, grp.indicator, scores = "Bartlett")
{
  latent.labels <- unique(grp.indicator) # Get unique latent variable labels

  # Create as matrix first, then coerce to data.frame. This is because it's easiest
  # to make a matrix with our desired dimensions, and then treat it as a data.frame
  # later on.
  latent.vars.mtx <- matrix(nrow = nrow(data), ncol = length(latent.labels))
  latent.vars.df <- as.data.frame(latent.vars.mtx)
  names(latent.vars.df) <- latent.labels

  for(label in latent.labels)
  {
    latent.cols.indic <- grp.indicator == label # TRUE if column corresponds to label
    latent.cols.df <- data[latent.cols.indic] # sub-data.frame of only those columns corresp. to label

    fa <- factanal(latent.cols.df, factors = 1, scores = scores) # perform factor analysis

    latent.vars.df[label] <- fa$scores[,1] # Get scores *as vector* so we don't update column name
  }

  return(latent.vars.df)
}

