#' A function to run a global test for significance of top chromosome results across a range of chromosome sizes
#'
#' This function runs a global test for significance of top chromosome results across a range of chromosome sizes.
#'
#' @param results.list A list of length d, where d is the number of chromosome sizes to be included in a global test.
#'  Each element of the list must itself be a list whose first element \code{observed.data} is a data.table containing
#'  the unique chromosome results from \code{combine.islands} for a given chromosome size. The second element \code{permutation.list}
#'  is a list containing all permutation results data.tables, again using the unique chromosome results output by \code{combine.islands}
#'  for each permutation.
#' @return A list containing the observed test statistic \code{obs.test.stat}, the p-value \code{pval},
#' and a vector of permutation test statistics \code{perm.test.stats}.
#' @importFrom matrixStats rowMaxs
#' @importFrom stats ecdf
#' @export

run.global.test <- function(results.list){

  #loop over chromosome sizes
  chrom.size.ks.list <- lapply(results.list, function(chrom.size.res){

    #error checking
    if (! "observed.data" %in% names(chrom.size.res)){

      stop("Each element of results.list must be a list containing an element titled observed.data.")
    }

    if (! "permutation.list" %in% names(chrom.size.res)){

      stop("Each element of results.list must be a list containing an element titled permutation.list.")
    }

    #grab the observed data
    obs.data <- chrom.size.res$observed.data
    obs.fitness.scores <- obs.data$fitness.score

    #compute the eCDF for the observed data
    obs.data.seq <- seq(min(obs.fitness.scores), max(obs.fitness.scores), length.out = 10000)
    obs.ecdf.fun <- ecdf(obs.fitness.scores)
    obs.ecdf <- obs.ecdf.fun(obs.data.seq)

    #compute the mean eCDF from the permutation results
    perm.list <- chrom.size.res$permutation.list
    perm.ecdf.mat <- matrix(NA, ncol = length(obs.data.seq), nrow = length(perm.list))
    for (i in 1:length(perm.list)){

      perm.res <- perm.list[[i]]
      perm.ecdf.fun <- ecdf(perm.res$fitness.score)
      perm.ecdf.mat[i, ] <- perm.ecdf.fun(obs.data.seq)

    }
    mean.perm.ecdf <- colMeans(perm.ecdf.mat)

    #now get the KS test stat compared with the mean eCDF
    obs.ks <- max(mean.perm.ecdf - obs.ecdf)
    mean.perm.ecdf.mat <- matrix(rep(mean.perm.ecdf, length(perm.list)), nrow = length(perm.list), byrow = TRUE)
    perm.ks <- rowMaxs(mean.perm.ecdf.mat - perm.ecdf.mat)
    return(list(obs.ks = obs.ks, perm.ks = perm.ks))

  })

  #grab the observed vector
  obs.ks.vec <- sapply(chrom.size.ks.list, function(x) x$obs.ks)

  #matrix of permutation based vectors
  perm.ks.mat <- sapply(chrom.size.ks.list, function(x) x$perm.ks)

  #mean permutation based vector
  mean.perm.vec <- colMeans(perm.ks.mat)

  #covariance for permutation based vector
  perm.cov.mat <- cov(perm.ks.mat)
  perm.cov.mat.inv <- solve(perm.cov.mat)

  #calculate mahalanobis distance to mean vector for obs data
  obs.mahala <- as.numeric((obs.ks.vec - mean.perm.vec) %*% perm.cov.mat.inv %*% (obs.ks.vec - mean.perm.vec))

  #now for each of the permutations
  perm.mahala <- rep(NA, nrow(perm.ks.mat))
  for (i in 1:nrow(perm.ks.mat)){

    perm.ks.vec <- perm.ks.mat[i, ]
    perm.mahala[i] <- (perm.ks.vec - mean.perm.vec) %*% perm.cov.mat.inv %*% (perm.ks.vec - mean.perm.vec)

  }

  #count the number of permutations greater than observed
  n.perms.greater <- sum(perm.mahala > obs.mahala)
  if (n.perms.greater == 0){

    pval <- 1/length(perm.mahala)
    warning(paste("No permutation Mahalanobis distances greater than observed, setting p-value to", pval))

  } else {

    pval <- n.perms.greater/length(perm.mahala)

  }

  #return results list
  res.list <- list(obs.test.stat = obs.mahala, pval = pval, perm.test.stats = perm.mahala)
  return(res.list)

}
