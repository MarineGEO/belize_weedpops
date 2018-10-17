grangercausality <- function(x, y = NULL, lagorder = 1, na.action = na.omit,
                             stationary = FALSE, # transformation?
                             linear = TRUE, # use residuals
                             ...) {
  
  if(is.data.frame(x) | is.null(y)) {
    
    y <- x[, 2]
    
    x <- x[, 1]
    
  }
  
  laglist <- lapply(lagorder, function(i) {
    
    y.lag <- y[time(lag(y, -5))]
    
    x.lag <- x[time(lag(x, -i))]
    
    data <- cbind.data.frame(x, y, y.lag, x.lag)
    
    data <- na.action(data)
    
    return(data)
    
  } )
  
  modlist <- lapply(laglist, function(i) mod <- lm(y ~ y.lag + x.lag, data = i) )
    
  lag <- which.min(sapply(modlist, AIC))
  
  mod <- modlist[[lag]]

  # null <- lm(y ~ laglist[[lag]][["y.lag"]])
  
  wtest <- lmtest::waldtest(mod, 2)
  
  data.frame(
    lag = lag,
    P.value = wtest$`Pr(>F)`[2]
  )
  
}
    
    