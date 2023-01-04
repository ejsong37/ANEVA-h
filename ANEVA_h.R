#' ANEVA-h
#'
#' Fits data to BLN distribution using Likelihood function then estimates vgs
#' 
#'
#' @param ref_table Dataframe containing reference count data
#' @param alt_table Dataframe containing alternate count data
#' @param min_occurances minimum number of people who contain the gene
#' @return Dataframe with fitted standard deviations, and coverage, and the cutoff to classify a sample as good quality.
#' 
ANEVA_h <- function(ref_table, alt_table, min_occurances = 10) {
  
  # Constants
  # Only analyze genes with expression higher than "MinExpression" in at least "(Mincases)" samples
  MINCASES <- 6 # Minimum number of samples which contain the gene
  MINEXPRESSION <- 30 # Minimum amount of expression of the gene
  STD_CHECKER <- 0 # Checks for BLN likelihood is defined at 0
  
  # Initializing Vectors
  gene_names <- c()
  vg_vec <- c()
  sigma_vec <- c()
  
  # Loops through all genes in table
  for(i in seq(ref_table$X)) { 
    # Gets all the counts into a specific row
    row_ref <- ref_table[i,]
    row_alt <- alt_table[i,]
    ref_counts <- c()
    alt_counts <- c()
    for(j in seq(3,length(ref_table))) { # Gets all ASE data for that one gene
      ref_counts <- c(ref_counts, row_ref[[j]])
      alt_counts <- c(alt_counts, row_alt[[j]])
    }
    
    total_counts <- ref_counts + alt_counts
    
    # Only analyze genes with expression higher than "MinExpression" in at least "(Mincases)" samples
    if(sum(total_counts > 0) >= min_occurances) {
      NThrSamples = sum(total_counts >= MINEXPRESSION) 
      if(NThrSamples >= MINCASES) {
        if(BLN_likelihood(ref_counts,total_counts, Std = STD_CHECKER) != Inf) { # Checks for BLN likelihood to be defined at 0
          if(sum(total_counts) > 5000) {
            sigma <- Fit_BLN(ref_counts, total_counts) # Optimizes for sigma (MLE)
            sigma_vec <- c(sigma_vec, sigma)
            vg <- estimate_vg(sigma) # Calculates vg for s
            # Saving to vectors
            gene_names <- c(gene_names, ref_table$name[i])
            vg_vec <- c(vg_vec, vg)
          }
        }
      }
    }
  }
  # Returns final result in list
  results <- list(GeneID = gene_names, vg = vg_vec)
  data_frame <- data.frame(results)
  return(data_frame)
}


