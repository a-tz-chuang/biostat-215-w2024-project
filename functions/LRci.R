#----- Function for Likelihood Ratio Based Confidence Intervals

# References:
# Thomas D.R., Grunkemeier, G.L. (1975) "Confidence Interval Estimation of Survival Probabilities for Censored Data".  Journal of the American Statistical Association, 70, 865 - 871.
# Li, G. (1995) "On nonparametric likelihood ratio estimation of survival probabilities for censored data". Statistics and Probability Letters, 25, 95 - 104.

LRci.surv <- function(fit, t, alpha = 0.05, upper.lim = 100, lower.lim = NULL) {
  
  #----- Inputs:
  # fit       : A survfit object
  # t         : Time at which we want to find confidence interval for
  # alpha     : (1 - alpha) * 100% confidence interval
  # upper.lim : Manually inputted value for upper limit of root (May need to be set larger if warning arises)
  # lower.lim : Manually inputted value for lower limit of root (default should be okay)
  
  rj <- fit$n.risk[fit$time <= t]  #Number at risk for tj < t
  dj <- fit$n.event[fit$time <= t] #NUmber of events for tj < t
  
  #- Function we want to find roots for
  f <- function(lambda) {
    -2 * sum((rj - dj) * log(1 + lambda / (rj - dj)) - rj * log(1 + lambda / rj)) -
               qchisq(1 - alpha, 1)
  } 
  
  if(is.null(lower.lim)) {
    lower.lim <- max(dj - rj) + 1
  }
  
  #- Find roots using uniroot
  lam_U <- uniroot(f = f, interval = c(0, upper.lim))
  lam_L <- uniroot(f = f, interval = c(lower.lim, 0))
  
  LL <- prod(1 - dj / (rj + lam_L$root))
  UL <- prod(1 - dj / (rj + lam_U$root))
  
  ci.lim <- cbind(LL, UL)
  
  #- List results here
  results <- list()
  
  results$time     <- t
  results$surv.est <- round(fit$surv[fit$time == max(fit$time[fit$time - t <= 0])], 3)
  results$conf.int <- round(ci.lim, 3)
  results$alpha    <- alpha
  return(results)
}
