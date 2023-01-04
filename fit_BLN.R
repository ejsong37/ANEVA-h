#' Fits Binomial Logit-Normal Model (BLN)
#'
#' @param ref_counts Reference read counts.
#' @param total_counts Total read counts.
#' @param start Optional Starting point for optimization (Mean, std, lambda)
#' @param factr relative tolerance
#' @param numCores Number of cores for parallelization
#' @return Estimated parameters generated from the fit

Fit_BLN <- function(ref_counts, total_counts) { 
  start.Std <- 0.5 # starting point for sigma prediction
  MLE_optim <- optim(par <- start.Std, fn = BLN_likelihood, ref_counts = ref_counts, total_counts = total_counts, lower = 0.01, upper = 1.75, method = "L-BFGS-B")  # Optimizer function
  return(MLE_optim$par)
}

#' Binomial Logit-Normal Likelihood Function
#'
#' @param ref_counts Reference read counts.
#' @param total_counts Total read counts.
#' @param std Standard Deviation
#' @return Likelihood
BLN_likelihood <- function(ref_counts,total_counts, Std) {
  L <- dbln(ref_counts, total_counts, 0, Std)
  NLL <- -sum(log(L)) # Negative log-likelihood 
  return(NLL)
}


