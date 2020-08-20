#' A function to run a global test of the null hypothesis that the observed fitness score CDF is not
#' more heavy-tailed on the right than the average null CDF across a range of chromosome sizes
#'
#' This function runs a global test across chromosome sizes of the null hypothesis
#' that the observed fitness score CDF is not more heavy-tailed on the right than the average null CDF
#'
#' @param results.list A list of length d, where d is the number of chromosome sizes to be included in a global test.
#'  Each element of the list must itself be a list whose first element \code{observed.data} is a vector of fitness scores from the
#'  the \code{all.results} chromosome results from \code{combine.islands} for a given chromosome size. The second element \code{permutation.list}
#'  is a list containing vectors of all permutation results fitness scores, again using the \code{all.results} results output by
#'  \code{combine.islands} for each permutation.
#' @param n.top.scores The number of top scoring chromosomes to be used in calculating the global test. Defaults to 1000.
#' @return A list containing the following:
#' \describe{
#'  \item{obs.test.stat}{The observed Mahalanobis distance global test statistic.}
#'  \item{pval}{The p-value for the global test.}
#'  \item{perm.test.stats}{A vector of Mahalanobis distance global test statistics for the permuted datasets.}
#'  \item{element.test.stats}{A vector of test statistics for the observed data, where each element of the vector
#'  corresponds to a specific chromosome size.}
#'  \item{element.pvals}{A vector of p-values corresponding to the test statistics in \code{element.test.stats}.}
#'  \item{perm.elem.test.stat.mat}{A matrix of test statistics for each of the permutation datasets.
#'  The rows of the matrix correspond to different permutations, and the columns correspond to chromosome sizes.}
#'  \item{obs.ks.vec}{A vector of observed Kolmogorov Smirnov test statistics for each chromosome size.}
#'  \item{perm.ks.mat}{A matrix of Kolmogorov Smirnov test statistics for the permutation datasets, where rows correspond to
#'  permutations and columns correspond to chromosome sizes.}
#'  \item{max.obs.fitness}{A vector of the maximum fitness score for each chromosome size in the observed data.}
#'  \item{max.perm.fitness}{A list of vectors for each chromosome size of maximum observed fitness scores for each permutation.}
#'  \item{max.order.pvals}{A vector of p-values for the maximum observed order statistics for each chromosome size.
#'   P-values are the proportion of permutation based maximum order statistics that exceed the observed maximum fitness score.}
#' }
#' @importFrom stats cov
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' data(snp.annotations)
#' library(Matrix)
#' set.seed(1400)
#' block.ld.mat <- as.matrix(bdiag(list(matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25))))
#'
# #preprocess data
#' pp.list <- preprocess.genetic.data(case[, 1:10], father.genetic.data = dad[ , 1:10],
#'                                mother.genetic.data = mom[ , 1:10],
#'                                block.ld.mat = block.ld.mat[ , 1:10])
#' ## run GA for observed data
#'
#' #observed data chromosome size 2
#' run.ga(pp.list, n.chromosomes = 5, chromosome.size = 2, results.dir = 'tmp_2',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  combined.res2 <- combine.islands('tmp_2', snp.annotations[ 1:10, ], pp.list)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#'  #observed data chromosome size 3
#'  run.ga(pp.list, n.chromosomes = 5, chromosome.size = 3, results.dir = 'tmp_3',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  combined.res3 <- combine.islands('tmp_3', snp.annotations[ 1:10, ], pp.list)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#' #create three permuted datasets
#' set.seed(1400)
#' perm.data.list <- permute.dataset(case[ , 1:10],
#'                                   father.genetic.data = dad[ , 1:10],
#'                                   mother.genetic.data = mom[ , 1:10],
#'                                   n.permutations = 3)
#'
#' #pre-process permuted data
#' p1.list <- preprocess.genetic.data(perm.data.list[['permutation1']]$case,
#'                                    complement.genetic.data = perm.data.list[['permutation1']]$comp,
#'                                    block.ld.mat = block.ld.mat[ , 1:10])
#'
#' p2.list <- preprocess.genetic.data(perm.data.list[['permutation2']]$case,
#'                                    complement.genetic.data = perm.data.list[['permutation2']]$comp,
#'                                    block.ld.mat = block.ld.mat[ , 1:10])
#'
#' p3.list <- preprocess.genetic.data(perm.data.list[['permutation3']]$case,
#'                                    complement.genetic.data = perm.data.list[['permutation3']]$comp,
#'                                    block.ld.mat = block.ld.mat[ , 1:10])
#'
#' ##run GA for permuted data
#'
#' #permutation 1, chromosome size 2
#' run.ga(p1.list, n.chromosomes = 5, chromosome.size = 2, results.dir = 'p1_tmp_2',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  p1.combined.res2 <- combine.islands('p1_tmp_2', snp.annotations[ 1:10, ], p1.list)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#' #permutation 1, chromosome size 3
#' run.ga(p1.list, n.chromosomes = 5, chromosome.size = 3, results.dir = 'p1_tmp_3',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  p1.combined.res3 <- combine.islands('p1_tmp_3', snp.annotations[ 1:10, ], p1.list)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#' #permutation 2, chromosome size 2
#' run.ga(p2.list, n.chromosomes = 5, chromosome.size = 2, results.dir = 'p2_tmp_2',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  p2.combined.res2 <- combine.islands('p2_tmp_2', snp.annotations[ 1:10, ], p2.list)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#' #permutation 2, chromosome size 3
#' run.ga(p2.list, n.chromosomes = 5, chromosome.size = 3, results.dir = 'p2_tmp_3',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  p2.combined.res3 <- combine.islands('p2_tmp_3', snp.annotations[ 1:10, ], p2.list)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#' #permutation 3, chromosome size 2
#' run.ga(p3.list, n.chromosomes = 5, chromosome.size = 2, results.dir = 'p3_tmp_2',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  p3.combined.res2 <- combine.islands('p3_tmp_2', snp.annotations[ 1:10, ], p3.list)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#' #permutation 3, chromosome size 3
#' run.ga(p3.list, n.chromosomes = 5, chromosome.size = 3, results.dir = 'p3_tmp_3',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  p3.combined.res3 <- combine.islands('p3_tmp_3', snp.annotations[ 1:10, ], p3.list)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#'  ## create list of results
#'
#'  #chromosome size 2 results
#'  chrom2.list <- list(observed.data = combined.res2$all.results$fitness.score,
#'                     permutation.list = list(p1.combined.res2$all.results$fitness.score,
#'                                             p2.combined.res2$all.results$fitness.score,
#'                                             p3.combined.res2$all.results$fitness.score))
#'
#'  #chromosome size 3 results
#'  chrom3.list <- list(observed.data = combined.res3$all.results$fitness.score,
#'                     permutation.list = list(p1.combined.res3$all.results$fitness.score,
#'                                             p2.combined.res3$all.results$fitness.score,
#'                                             p3.combined.res3$all.results$fitness.score))
#'
#'  final.results <- list(chrom2.list, chrom3.list)
#'
#'  ## run global test
#'  global.test.res <- run.global.test(final.results)
#'
#'  lapply(c('tmp_2', 'tmp_3', 'p1_tmp_2', 'p2_tmp_2', 'p3_tmp_2',
#'           'p1_tmp_3', 'p2_tmp_3', 'p3_tmp_3'), unlink, recursive = TRUE)
#'
#'
#' @importFrom matrixStats rowMaxs
#' @importFrom stats ecdf
#' @export

run.global.test <- function(results.list, n.top.scores = 1000) {

    # loop over chromosome sizes
    chrom.size.ks.list <- lapply(results.list, function(chrom.size.res) {

        # error checking
        if (!"observed.data" %in% names(chrom.size.res)) {

            stop("Each element of results.list must be a list containing an element titled observed.data.")
        }

        if (!"permutation.list" %in% names(chrom.size.res)) {

            stop("Each element of results.list must be a list containing an element titled permutation.list.")
        }

        # grab the observed data
        obs.fitness.scores <- sort(chrom.size.res$observed.data, decreasing = TRUE)[seq_len(n.top.scores)]

        # grab permuted data
        perm.list <- lapply(chrom.size.res$permutation.list, function(perm.scores){

            sort(perm.scores, decreasing = TRUE)[seq_len(n.top.scores)]

        })

        # get list of all unique fitness scores
        all.fitness.scores <- sort(unique(c(obs.fitness.scores, unlist(perm.list))))

        # compute the eCDF for the observed data
        obs.ecdf.fun <- ecdf(obs.fitness.scores)
        obs.ecdf <- obs.ecdf.fun(all.fitness.scores)

        # compute ecdf for each permutation
        perm.ecdf.list <- lapply(perm.list, function(permutation){

            ecdf(permutation)

        })

        # compute the mean eCDF from the permutation results
        perm.ecdf.mat <- matrix(NA, ncol = length(all.fitness.scores), nrow = length(perm.list))
        for (i in seq_len(length(perm.list))) {

            perm.ecdf.fun <- perm.ecdf.list[[i]]
            perm.ecdf.mat[i, ] <- perm.ecdf.fun(all.fitness.scores)

        }
        mean.perm.ecdf <- colMeans(perm.ecdf.mat)

        # now get the KS test stat compared with the mean eCDF
        obs.ks <- max(mean.perm.ecdf - obs.ecdf)

        # now for each of the permutations
        mean.perm.ecdf.mat <- matrix(rep(mean.perm.ecdf, length(perm.list)),
                                     nrow = length(perm.list), byrow = TRUE)
        perm.ks <- rowMaxs(mean.perm.ecdf.mat - perm.ecdf.mat)

        return(list(obs.ks = obs.ks, perm.ks = perm.ks))

    })

    # grab the observed vector
    obs.ks.vec <- vapply(chrom.size.ks.list, function(x) x$obs.ks, 1.0)

    # matrix of permutation based vectors
    perm.ks.mat <- vapply(chrom.size.ks.list, function(x) x$perm.ks,
                          rep(1.0, length(chrom.size.ks.list[[1]]$perm.ks)))

    # mean permutation based vector
    mean.perm.vec <- colMeans(perm.ks.mat)

    # covariance for permutation based vector
    perm.cov.mat <- cov(perm.ks.mat)
    perm.cov.mat.inv <- solve(perm.cov.mat)

    # calculate mahalanobis distance to mean vector for obs data
    obs.mean.dif <- obs.ks.vec - mean.perm.vec
    obs.mean.dif[obs.mean.dif < 0] <- 0
    obs.mahala <- as.numeric(obs.mean.dif %*% perm.cov.mat.inv %*% obs.mean.dif)
    #obs.mahala <- sqrt(sum(obs.mean.dif^2))
    #obs.mahala <- sqrt(sum((obs.mean.dif*diag(perm.cov.mat.inv))^2))

    # now for each of the permutations
    perm.mahala <- rep(NA, nrow(perm.ks.mat))
    for (i in seq_len(nrow(perm.ks.mat))) {

        perm.ks.vec <- perm.ks.mat[i, ]
        perm.mean.dif <- perm.ks.vec - mean.perm.vec
        perm.mean.dif[perm.mean.dif < 0] <- 0
        #perm.mahala[i] <- sqrt(sum(perm.mean.dif^2))
        perm.mahala[i] <- perm.mean.dif %*% perm.cov.mat.inv %*% perm.mean.dif
        #perm.mahala[i] <- sqrt(sum((perm.mean.dif*diag(perm.cov.mat.inv))^2))
    }

    # count the number of permutations greater than observed
    n.perms.greater <- sum(perm.mahala > obs.mahala)
    if (n.perms.greater == 0) {

        pval <- 0
        suggested <- paste0("1/", length(perm.mahala))
        print(paste("No permutation Mahalanobis distances greater than observed, consider setting p-value to",
            suggested))

    } else {

        pval <- n.perms.greater/length(perm.mahala)

    }

    # also look at element-wise results
    obs.elements <- as.numeric((obs.ks.vec - mean.perm.vec)/diag(perm.cov.mat.inv))^2
    perm.elem.mat <- matrix(NA, nrow(perm.ks.mat), ncol = ncol(perm.ks.mat))
    for (i in seq_len(nrow(perm.ks.mat))) {

        perm.ks.vec <- perm.ks.mat[i, ]
        perm.elem.mat[i, ] <- ((perm.ks.vec - mean.perm.vec)/diag(perm.cov.mat.inv))^2

    }

    # element wise p-vals
    element.pvals <- vapply(seq_len(length(obs.elements)), function(element.position) {

        obs.elem <- obs.elements[element.position]
        perm.elems <- perm.elem.mat[, element.position]
        pval <- sum(perm.elems > obs.elem)/length(perm.elems)
        return(pval)

    }, 1.0)

    # maximum permutation based fitness scores
    max.perm.fitness <- vapply(results.list, function(chrom.size.res) {

        perm.list <- chrom.size.res$permutation.list
        vapply(perm.list, max, 1.0)

    }, rep(1.0, length(results.list[[1]]$permutation.list)))

    # maximum observed fitness scores
    max.obs.fitness <- vapply(results.list, function(chrom.size.res) {

        obs.fitness.scores <- chrom.size.res$observed.data
        max(obs.fitness.scores)

    }, 1.0)

    # pvals for max order statistics
    max.order.pvals <- vapply(seq_along(max.obs.fitness), function(chrom.size) {

        max.obs <- max.obs.fitness[chrom.size]
        max.perms <- max.perm.fitness[ , chrom.size]
        pval <- sum(max.perms > max.obs)/length(max.perms)
        pval

    }, 1.0)

    # return results list
    res.list <- list(obs.test.stat = obs.mahala, pval = pval, perm.test.stats = perm.mahala, element.test.stats = obs.elements,
        element.pvals = element.pvals, perm.elem.test.stat.mat = perm.elem.mat, obs.ks.vec = obs.ks.vec,
        perm.ks.mat = perm.ks.mat, max.obs.fitness = max.obs.fitness, max.perm.fitness = max.perm.fitness,
        max.order.pvals = max.order.pvals)
    return(res.list)

}
