cbands.region <- function(surv.object, tL, tU) {
  
  #Error checking
  if(tL < min(surv.object$time)) stop("Lower limit is smaller than smallest observed survival time, choose larger value.")
  n <- length(surv.object$time)
  lower.lim <- max(surv.object$time[surv.object$time <= tL])
  upper.lim <- max(surv.object$time[surv.object$time <= tU])
  
  aL <- (n * surv.object$std.err[surv.object$time == lower.lim]^2) / (1 + n * surv.object$std.err[surv.object$time == lower.lim]^2)
  aU <- (n * surv.object$std.err[surv.object$time == upper.lim]^2) / (1 + n * surv.object$std.err[surv.object$time == upper.lim]^2)
  
  writeLines("Find critical regions in Klein and Moeschberger 2nd ed. (Appendix C.3a - C.4c)")
  results <- list()
  results$aL <- round(aL, 1)
  results$aU <- round(aU, 1)
  return(results)
}   

cbands.interval <- function(surv.object, tL, tU, crit.value, alpha = 0.05, type = "linear", method = "ep") {
  if(type != "linear" & type != "log" & type != "asin") stop("type must be one of the three options: linear, log, asin")
  if(method != "ep" & method != "hw") stop("method must be either 'ep' or 'hw'")  
  new.time <- surv.object$time[(surv.object$time >= max(surv.object$time[surv.object$time <= tL])) 
                               & (surv.object$time <= max(surv.object$time[surv.object$time <= tU]))]
  
  new.surv <-  surv.object$surv[(surv.object$time >= max(surv.object$time[surv.object$time <= tL])) 
                                & (surv.object$time <= max(surv.object$time[surv.object$time <= tU]))]
  
  new.se <- surv.object$std.err[(surv.object$time >= max(surv.object$time[surv.object$time <= tL])) 
                                & (surv.object$time <= max(surv.object$time[surv.object$time <= tU]))]
  
  new.event <- surv.object$n.event[(surv.object$time >= max(surv.object$time[surv.object$time <= tL])) 
                                   & (surv.object$time <= max(surv.object$time[surv.object$time <= tU]))]
  
  results <- data.frame(t = new.time, surv = new.surv, se = new.se)  
  #- Equal Precision Bands
  if(method == "ep") {
    if(type == "linear") {
      results$LL <- new.surv - crit.value * new.se * new.surv
      results$UL <- new.surv + crit.value * new.se * new.surv
    } else if (type == "log") {
      theta <- exp(crit.value * new.se / log(new.surv))
      results$LL <- new.surv^(1 / theta)
      results$UL <- new.surv^theta
    } else if (type == "asin") {
      results$LL <- sin(sapply(asin(sqrt(new.surv)), function(x) max(0, x)) - 0.5 * crit.value * new.se * 
                          sqrt(new.surv / (1 - new.surv)))^2
      results$UL <- sin(sapply(asin(sqrt(new.surv)), function(x) min(pi / 2, x)) + 0.5 * crit.value * new.se * 
                          sqrt(new.surv / (1 - new.surv)))^2 
    }
  } else if (method == "hw") { #Hall-Wellner Bands
    n <- length(surv.object$time)
    if(type == "linear") {
      results$LL <- new.surv - crit.value * (1 + n * new.se^2) / sqrt(n) * new.surv
      results$UL <- new.surv + crit.value * (1 + n * new.se^2) / sqrt(n) * new.surv
    } else if (type == "log") {
      theta <- exp(crit.value * (1 + n * new.se^2) / (sqrt(n) * log(new.surv)))
      results$LL <- new.surv^(1 / theta)
      results$UL <- new.surv^theta
    } else if (type == "asin") {
      results$LL <- sin(sapply(asin(sqrt(new.surv)), function(x) max(0, x)) - 0.5 * crit.value * (1 + n * new.se^2) / sqrt(n) * 
                          sqrt(new.surv / (1 - new.surv)))^2
      results$UL <- sin(sapply(asin(sqrt(new.surv)), function(x) min(pi / 2, x)) + 0.5 * crit.value * (1 + n * new.se^2) / sqrt(n) * 
                          sqrt(new.surv / (1 - new.surv)))^2 
    }
  }
  
  #----- Returning results here
  writeLines(paste0("Returning ", type, "-type confidence bands using ", method, " method."))
  return(round(results[which(new.event != 0), ], 3))
}

#cbands.interval(surv.object, tL = 100, tU = 600, crit.value = 1.3211, type = "linear", method = "hw")
