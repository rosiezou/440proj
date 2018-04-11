library(truncnorm)
library(MASS)

# Generate imputed response
gen.imp.resp <- function(df, num.iter = 20)
{
  df.interp <- df # For initial interpolation
  
  for(q in names(df))
  {
    responses <- df[[q]]
    min.resp <- min(responses)
    max.resp <- max(responses)
    
    probs <- cumsum(colSums(sapply(min.resp:max.resp, 
                                   FUN = function(i) {responses == i})) / nrow(df))
    
    quantiles <- c(-Inf, qnorm(probs)) # Get all quantiles at the probs, including -Inf
    
    lower <- quantiles[responses] # Get lower bound of bins for responses
    upper <- quantiles[responses + 1] # Same but get upper bound
    
    # Generate from truncated normal
    sim.truncnorm <- rtruncnorm(n = 1, a = lower, b = upper, mean = 0, sd = 1)
    
    df.interp[q] <- sim.truncnorm
  }
  
  df.analysis <- df.interp
  
  for (i in 1:num.iter) # number of updates
  {
    for(q in names(df.analysis)) # For each column
    {
      df.analysis[q] <- df[[q]]
      
      # Construct regression formula with q as the response
      regression.formula <- as.formula(paste("as.ordered(", q, ")", "~ ."))
      inhs.polr <- polr(regression.formula, data = df.analysis, method = "probit")
      quantiles <- c(-Inf, inhs.polr$zeta, Inf) # Get quantiles that include -Inf and Inf
      
      sim.bins <- as.numeric(simulate(inhs.polr)[[1]]) # Simulate new bins
      
      lower <- quantiles[sim.bins] # Get lower bounds of bins in sim.bins
      upper <- quantiles[sim.bins + 1] # Get upper bounds of bins in sim.bins
      
      # Generate from truncated normal
      sim.truncnorm <- rtruncnorm(n = 1, a = lower, b = upper, mean = 0, sd = 1)
      
      df.analysis[q] <- sim.truncnorm # Update column
    } 
  }
  
  return(df.analysis)
}

# Test 
setwd("~/GitProjects/440proj") # Set to project folder
multiis <- read.csv("data/MULTIIS.csv")
imputed <- gen.imp.resp(multiis)

head(imputed)

# df is the data frame of on which to perform factor analysis
# grp.indic is the vector of length (# columns in df) specifying which latent variable
# each column refers to.
gen.latent.vars <- function(df, grp.indicator, scores = "Bartlett")
{
  latent.labels <- unique(grp.indicator) # Get unique latent variable labels
  
  # Create as matrix first, then coerce to data.frame. This is because it's easiest
  # to make a matrix with our desired dimensions, and then treat it as a data.frame
  # later on. 
  latent.vars.mtx <- matrix(nrow = nrow(df), ncol = length(latent.labels))
  latent.vars.df <- as.data.frame(latent.vars.mtx) 
  names(latent.vars.df) <- latent.labels
  
  for(label in latent.labels)
  {
    latent.cols.indic <- grp.indicator == label # TRUE if column corresponds to label
    latent.cols.df <- df[latent.cols.indic] # sub-data.frame of only those columns corresp. to label
    
    fa <- factanal(latent.cols.df, factors = 1, scores = scores) # perform factor analysis
    
    latent.vars.df[label] <- fa$scores[,1] # Get scores *as vector* so we don't update column name
  }
  
  return(latent.vars.df)
}

# Test
# Create indicators (a label indicating which latent variable the question corresponds to)
grp.indicator <- sapply(names(imputed), FUN = 
                          function(x){strsplit(x, split = "_")[[1]][2]})

latent.vars <- gen.latent.vars(imputed, grp.indicator = grp.indicator)

head(latent.vars)
  
  