#' Estimates the expected variance in dosage that is due to inter-individual
#' genetic differences in population (VG)
#'
#' @param std standard deviation
#' @return VG estimate
estimate_vg <- function(std) {
  # This function calculates the VG estimates of the dataset
  # Parameters: sigma, sr, data
  # sigma is the standard deviation estimated from earlier optimization
  # sr is the residual regulatory variance
  vg_estimates = integrate(formula_vg(std), lower = -4*std,upper = 4*std,std = std)$value
  return(vg_estimates)
}

#' Helper function for estimate_vg
#'
#' @param sigma Reference read counts.
#' @return Estimated parameters generated from the fit
formula_vg <- function(sigma) {
  # This function gives the formula for the calculation of the VG estimates to be integrated
  # Parameters: sigma, sr
  # sigma is the standard deviation estimated from earlier optimization
  # sr is the residual regulatory variance
  inner_formula <- function(sigma,sr) {
    expect <- formula_expectation(sigma) #formula for formula expectation
    expectation = integrate(expect,-4*sigma,4*sigma,sigma = sigma)$value
    formula = ((log(exp(sr)+1) - expectation)**2)*p_sr(sigma,sr)
    return(formula)
  }
  return(inner_formula)
}

#' Helper function for formula_vg
#'
#' @param ref_counts Reference read counts.
#' @param total_counts Total read counts.
#' @param start Optional Starting point for optimization (Mean, std, lambda)
#' @param factr relative tolerance
#' @param numCores Number of cores for parallelization
#' @return Estimated parameters generated from the fit
formula_expectation <- function(sigma) {
  # This function gives the expectation of ln(e) to be used in formula_vg
  # Parameters: sigma, sr
  # sigma is the standard deviation estimated from earlier optimization
  # sr is the residual regulatory variance
  inner_formula <- function(sigma,sr) {
    formula = log(exp(sr)+1)*p_sr(sigma,sr)
    return(formula)
  }
  return(inner_formula)
}

#' Calculates the probability of sr used in formula_vg
#'
#' @param std Standard Deviation
#' @param sr Residual Regulatory variance
#' @return Estimated parameters generated from the fit
p_sr <- function(sigma, sr) {
  return(dnorm(sr, 0, sigma)*(1/0.999))
}

