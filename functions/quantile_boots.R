quantile_boots <- function(t1, c1, q = .5, R = 1000, seed = 1991, plots = FALSE){
  
  #----- Inputs:
  # t1    : Vector of observed survival times
  # c1    : Vector of censoring indicators (0 = censored).
  # q     : Quantile of interest to be compared (Default = .5).
  # R     : Number of iterations for bootstrapping. [Higher = More Time]
  # seed  : Setting a seed for reproducibility. 
  # plots : Whether or not to plot the KM curve
  
  prob     <- 1 - q
  fit      <- survfit(Surv(t1,c1) ~ 1, conf.type = 'none') # Control
  q.est    <- unname(quantile(fit, prob = q))
  if(is.na(q.est)){
    stop(paste("Estimated survival time for ", q*100, "th quantile was not found. Program Stopped."))
  }
  
  Qp <- function(t1, c1){
    fit1 <- survfit(Surv(t1,c1) ~ 1, conf.type = 'none') # Control
    Finv <- unname(quantile(fit1, prob = q)) #Finv(p) 
    if(is.na(Finv)){
      warning("Median survival time could not be estimated for sample. Largest observed survival time was used as median.")
      Finv <- max(t1)
    }
    return(Finv)
  }
  
  bootstrap.est <- numeric(R)
  for(i in 1:R){
    btsp1 <- sample(c(1:length(t1)), replace = TRUE)
    t1.bootstrap <- t1[btsp1]
    c1.bootstrap <- c1[btsp1]
    bootstrap.est[i] <- Qp(t1.bootstrap, c1.bootstrap) #Bootstrapped KM
  }
  
  if(plots == TRUE){
    plot(fit, col = 'red', ylab = 'Estimated Survival Function', xlab = 'Time', main = 'Kaplan-Meier Estimate')
    lines(x=c(q.est, q.est), y = c(-.5, prob), lty = 2)
    lines(x=c(0, q.est), y = c(prob, prob), lty = 2)
  }
  
  return(list("Quantile Estimate" = q.est, "Std. Err" = sd(bootstrap.est)))
  
  #----- Output
  # Quantile Estimate: Estimate of the quantile using KM.
  # Std. Err: Standard error of quantle estimate from bootstrap.
  # A Kaplan-Meier curve is also produced (if wanted).
}
