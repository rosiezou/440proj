# Generate imputed response
gen.imp.resp <- function(df)
{
  df.interp <- df # For initial interpolation
  
  for(q in names(df))
  {
    responses <- df[[q]]
    min.resp <- min(responses)
    max.resp <- max(responses)
    
    probs <- cumsum(colSums(sapply(min.resp:max.resp, 
                                   FUN = function(i) {responses == i})) / nrow(df))
    
    quantiles <- c(-Inf, qnorm(probs))
    
    lower <- quantiles[responses] # Get lower bound of bins for responses
    upper <- quantiles[responses + 1] # Same but get upper bound
    
    # Generate from truncated normal
    sim.truncnorm <- rtruncnorm(n = 1, a = lower, b = upper, mean = 0, sd = 1)
    
    df.interp[q] <- sim.truncnorm
  }
  
  df.analysis <- df.interp
  
  for(q in names(df.analysis)) # For each column
  {
    df.analysis[q] <- df[[q]]
    
    # Construct regression formula with q as the response
    regression.formula <- as.formula(paste("as.ordered(", q, ")", "~ ."))
    inhs.polr <- polr(regression.formula, data = df.analysis, method = "probit")
    
    sim.bins <- as.numeric(simulate(inhs.polr)[[1]]) # Simulate new bins
    
    lower <- zeta[sim.bins] # Get lower bounds of bins in sim.bins
    upper <- zeta[sim.bins + 1] # Get upper bounds of bins in sim.bins
    
    # Generate from truncated normal
    sim.truncnorm <- rtruncnorm(n = 1, a = lower, b = upper, mean = 0, sd = 1)
    
    df.analysis[q] <- sim.truncnorm # Update column
  }
  
  return(df.analysis)
}

# Test 
setwd("~/GitProjects/440proj") # Set to project folder
multiis <- read.csv("data/MULTIIS.csv")
imputed <- gen.imp.resp(multiis)

head(imputed)
