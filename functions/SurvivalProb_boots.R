SurvivalProb_boots <- function(t1, c1, t, R = 1000, seed = 1991, plots = FALSE) {
  
  fit <- survfit(Surv(t1,c1) ~ 1, conf.type = 'none') # Control
  if (t > max(t1) | t < min(t1)) {
    stop("The input time is beyond the time range of samples. Please try again!")
  }
  
  a <- summary(fit, times = t)
  point.est <- a$surv
  if (is.na(point.est)) {
    stop("Something wrong with the estimation or input. Please try again!")
  }
  
  Qp <- function(t1, c1){
    fit1 <- survfit(Surv(t1,c1) ~ 1, conf.type = 'none') # Control
    a <- summary(fit1, times = t)
    aa <- a$surv
    if (is.na(aa)) {
      stop("Something wrong with the estimation or input. Please try again!")
    } else {
      return(aa)
    }
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
    lines(x=c(t, t), y = c(-.5, point.est), lty = 2)
    lines(x=c(0, t), y = c(point.est, point.est), lty = 2)
  }
  
  return(list("Time" = t, "Survival Probability Estimate" = point.est, "Std. Err" = sd(bootstrap.est)))
}