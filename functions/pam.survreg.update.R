pam.survreg <- function(fit.survreg) 
{
  x.matrix.unsorted <- fit.survreg$x
  y.unsorted <- fit.survreg$y[, 1]
  censor.unsorted <- fit.survreg$y[, 2]
  nsize <- length(y.unsorted)
  y <- sort(y.unsorted)
  delta <- censor.unsorted[order(y.unsorted)]
  p <- dim(as.matrix(x.matrix.unsorted))[2]
  if (p == 1) {
    x.matrix <- as.matrix(x.matrix.unsorted[order(y.unsorted)])
  }
  else {
    x.matrix <- x.matrix.unsorted[order(y.unsorted), ]
  }
  nsize <- length(y)
  fit.km.censoring <- survfit(Surv(y, 1 - delta) ~ 1)
  sum.km.censoring <- summary(fit.km.censoring, times = y, 
                              extend = TRUE)
  km.censoring <- sum.km.censoring$surv
  km.censoring.minus <- c(1, km.censoring[-length(km.censoring)])
  ratio.km <- delta/km.censoring.minus
  ratio.km[is.nan(ratio.km)] <- 0
  weight.km <- ratio.km/(sum(ratio.km))
  if (fit.survreg$dist == "exponential") {
    t.predicted <- predict(fit.survreg, newdata = data.frame(x.matrix), 
                           type = "response") * gamma(2)
  }
  else if (fit.survreg$dist == "weibull") {
    t.predicted <- predict(fit.survreg, newdata = data.frame(x.matrix), 
                           type = "response") * gamma(1 + fit.survreg$scale)
  }
  else if (fit.survreg$dist == "lognormal") {
    t.predicted <- predict(fit.survreg, newdata = data.frame(x.matrix), 
                           type = "response") * exp((fit.survreg$scale)^2/2)
  }
  else if (fit.survreg$dist == "loglogistic") {
    t.predicted <- predict(fit.survreg, newdata = data.frame(x.matrix), 
                           type = "response") * gamma(1 + fit.survreg$scale) * 
      gamma(1 - fit.survreg$scale)
  }
  else {
    t.predicted <- predict(fit.survreg, newdata = data.frame(x.matrix), 
                           type = "response")
  }
  wls.fitted <- tryCatch(lm(y ~ t.predicted, weights = weight.km), 
                         error = function(e) {
                           return(c(NA, NA))
                         })
  calibrate.fitted <- tryCatch(predict(wls.fitted), error = function(e) {
    return(c(NA, NA))
  })
  num.rho2 <- sum(weight.km * (calibrate.fitted - sum(weight.km * 
                                                        y))^2)
  denom.rho2 <- sum(weight.km * (y - sum(weight.km * y))^2)
  R2 <- format(round(num.rho2/denom.rho2, digits = 4), nsmall = 4)
  num.L2 <- sum(weight.km * (y - calibrate.fitted)^2)
  denom.L2 <- sum(weight.km * (y - t.predicted)^2)
  L2 <- format(round(num.L2/denom.L2, digits = 4), nsmall = 4)
  return(list(R.squared = R2, L.squared = L2))
}
