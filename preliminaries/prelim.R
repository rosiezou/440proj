# Generate imputed response
gen.imp.resp <- function(df)
{
  df.analysis <- df
  
  for(q in names(df.analysis))
  {
    regression.formula <- as.formula(paste("as.ordered(", q, ")", "~ ."))
    inhs.polr <- polr(regression.formula, data = df.analysis, method = "probit")
    
    sim.bins <- as.numeric(simulate(inhs.polr)[[1]])
    
    lower <- zeta[sim.bins]
    upper <- zeta[sim.bins + 1]
    
    sim.truncnorm <- rtruncnorm(n = 1, a = lower, b = upper, mean = 0, sd = 1)
    
    df.analysis[q] <- sim.truncnorm
  }
  
  df.analysis
}