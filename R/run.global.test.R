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
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' library(Matrix)
#' chrom.mat <- as.matrix(bdiag(list(matrix(rep(TRUE, 2500^2), nrow = 2500),
#'                               matrix(rep(TRUE, 2500^2), nrow = 2500),
#'                               matrix(rep(TRUE, 2500^2), nrow = 2500),
#'                               matrix(rep(TRUE, 2500^2), nrow = 2500))))
#'
#  #preprocess data
#' pp.list <- preprocess.genetic.data(case[, 1:10], father.genetic.data = dad[ , 1:10],
#'                                mother.genetic.data = mom[ , 1:10],
#'                                chrom.mat = chrom.mat[ , 1:10])
#' ## run GA for observed data
#'
#' #observed data chromosome size 2
#' run.ga(pp.list, n.chromosomes = 5, chromosome.size = 2, results.dir = "tmp_2",
#'        cluster.type = "interactive", registryargs = list(file.dir = "tmp_reg", seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  combined.res2 <- combine.islands("tmp_2")
#'  unlink("tmp_reg", recursive = TRUE)
#'
#'  #observed data chromosome size 3
#'  run.ga(pp.list, n.chromosomes = 5, chromosome.size = 3, results.dir = "tmp_3",
#'        cluster.type = "interactive", registryargs = list(file.dir = "tmp_reg", seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  combined.res3 <- combine.islands("tmp_3")
#'  unlink("tmp_reg", recursive = TRUE)
#'
#' #create three permuted datasets
#' set.seed(1400)
#' perm.data.list <- permute.dataset(case[ , 1:10],
#'                                   father.genetic.data = dad[ , 1:10],
#'                                   mother.genetic.data = mom[ , 1:10],
#'                                   n.permutations = 3)
#'
#' #pre-process permuted data
#' p1.list <- preprocess.genetic.data(perm.data.list[["permutation1"]]$case,
#'                                    complement.genetic.data = perm.data.list[["permutation1"]]$comp,
#'                                    chrom.mat = chrom.mat[ , 1:10])
#'
#' p2.list <- preprocess.genetic.data(perm.data.list[["permutation2"]]$case,
#'                                    complement.genetic.data = perm.data.list[["permutation2"]]$comp,
#'                                    chrom.mat = chrom.mat[ , 1:10])
#'
#' p3.list <- preprocess.genetic.data(perm.data.list[["permutation3"]]$case,
#'                                    complement.genetic.data = perm.data.list[["permutation3"]]$comp,
#'                                    chrom.mat = chrom.mat[ , 1:10])
#'
#' ##run GA for permuted data
#'
#' #permutation 1, chromosome size 2
#' run.ga(p1.list, n.chromosomes = 5, chromosome.size = 2, results.dir = "p1_tmp_2",
#'        cluster.type = "interactive", registryargs = list(file.dir = "tmp_reg", seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  p1.combined.res2 <- combine.islands("p1_tmp_2")
#'  unlink("tmp_reg", recursive = TRUE)
#'
#' #permutation 1, chromosome size 3
#' run.ga(p1.list, n.chromosomes = 5, chromosome.size = 3, results.dir = "p1_tmp_3",
#'        cluster.type = "interactive", registryargs = list(file.dir = "tmp_reg", seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  p1.combined.res3 <- combine.islands("p1_tmp_3")
#'  unlink("tmp_reg", recursive = TRUE)
#'
#' #permutation 2, chromosome size 2
#' run.ga(p2.list, n.chromosomes = 5, chromosome.size = 2, results.dir = "p2_tmp_2",
#'        cluster.type = "interactive", registryargs = list(file.dir = "tmp_reg", seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  p2.combined.res2 <- combine.islands("p2_tmp_2")
#'  unlink("tmp_reg", recursive = TRUE)
#'
#' #permutation 2, chromosome size 3
#' run.ga(p2.list, n.chromosomes = 5, chromosome.size = 3, results.dir = "p2_tmp_3",
#'        cluster.type = "interactive", registryargs = list(file.dir = "tmp_reg", seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  p2.combined.res3 <- combine.islands("p2_tmp_3")
#'  unlink("tmp_reg", recursive = TRUE)
#'
#' #permutation 3, chromosome size 2
#' run.ga(p3.list, n.chromosomes = 5, chromosome.size = 2, results.dir = "p3_tmp_2",
#'        cluster.type = "interactive", registryargs = list(file.dir = "tmp_reg", seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  p3.combined.res2 <- combine.islands("p3_tmp_2")
#'  unlink("tmp_reg", recursive = TRUE)
#'
#' #permutation 3, chromosome size 3
#' run.ga(p3.list, n.chromosomes = 5, chromosome.size = 3, results.dir = "p3_tmp_3",
#'        cluster.type = "interactive", registryargs = list(file.dir = "tmp_reg", seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  p3.combined.res3 <- combine.islands("p3_tmp_3")
#'  unlink("tmp_reg", recursive = TRUE)
#'
#'  ## create list of results
#'
#'  #chromosome size 2 results
#'  chrom2.list <- list(observed.data = combined.res2$unique.results,
#'                     permutation.list = list(p1.combined.res2$unique.results,
#'                                             p2.combined.res2$unique.results,
#'                                             p3.combined.res2$unique.results))
#'
#'  #chromosome size 3 results
#'  chrom3.list <- list(observed.data = combined.res3$unique.results,
#'                     permutation.list = list(p1.combined.res3$unique.results,
#'                                             p2.combined.res3$unique.results,
#'                                             p3.combined.res3$unique.results))
#'
#'  final.results <- list(chrom2.list, chrom3.list)
#'
#'  ## run global test
#'  global.test.res <- run.global.test(final.results)
#'
#'  lapply(c("tmp_2", "tmp_3", "p1_tmp_2", "p2_tmp_2", "p3_tmp_2",
#'           "p1_tmp_3", "p2_tmp_3", "p3_tmp_3"), unlink, recursive = TRUE)
#'
#'
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
    obs.data.seq <- seq(0, max(obs.fitness.scores), length.out = 10000)
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

    pval <- 0
    suggested <- 1/length(perm.mahala)
    print(paste("No permutation Mahalanobis distances greater than observed, consider setting p-value to", suggested))

  } else {

    pval <- n.perms.greater/length(perm.mahala)

  }

  #also look at element-wise results
  obs.elements <- as.numeric((obs.ks.vec - mean.perm.vec) / diag(perm.cov.mat.inv))^2
  perm.elem.mat <- matrix(NA, nrow(perm.ks.mat), ncol = ncol(perm.ks.mat))
  for (i in 1:nrow(perm.ks.mat)){

    perm.ks.vec <- perm.ks.mat[i, ]
    perm.elem.mat[i, ] <- ((perm.ks.vec - mean.perm.vec) / diag(perm.cov.mat.inv))^2

  }

  #element wise p-vals
  element.pvals <- sapply(1:length(obs.elements), function(element.position){

    obs.elem <- obs.elements[element.position]
    perm.elems <- perm.elem.mat[ , element.position]
    pval <- sum(perm.elems > obs.elem)/length(perm.elems)
    return(pval)

  })

  #return results list
  res.list <- list(obs.test.stat = obs.mahala, pval = pval, perm.test.stats = perm.mahala,
                   element.test.stats = obs.elements, element.pvals = element.pvals,
                   perm.elem.test.stat.mat = perm.elem.mat)
  return(res.list)

}
