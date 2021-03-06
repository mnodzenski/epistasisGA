#' A function to run a global test of the null hypothesis that there are no joint SNP-disease associations
#' across a range of chromosome sizes
#'
#' This function runs a global test of the null hypothesis that there are no joint SNP-disease associations
#' across a range of chromosome sizes
#'
#' @param results.list A list of length d, where d is the number of chromosome sizes to be included in a global test.
#'  Each element of the list must itself be a list whose first element \code{observed.data} is a vector of fitness scores
#'  from \code{combine.islands} for a given chromosome size. The second element \code{permutation.list}
#'  is a list containing vectors of all permutation results fitness scores, again using the results output by
#'  \code{combine.islands} for each permutation.
#' @param n.top.scores The number of top scoring chromosomes, for each chromosome size, to be used in calculating the global test. Defaults to 30.
#' @return A list containing the following:
#' \describe{
#'  \item{obs.test.stat}{The observed test statistic.}
#'  \item{perm.test.stats}{A vector of test statistics from permuted data.}
#'  \item{pval}{The p-value for the global test.}
#'  \item{obs.marginal.test.stats}{A vector of observed test statistics for each chromosome size.}
#'  \item{perm.marginal.test.stats.mat}{A matrix of test statistics for the permutation datasets, where rows correspond to
#'  permutations and columns correspond to chromosome sizes.}
#'  \item{marginal.pvals}{A vector containing marignal p-values for each chromosome size.}
#'  \item{max.obs.fitness}{A vector of the maximum fitness score for each chromosome size in the observed data.}
#'  \item{max.perm.fitness}{A list of vectors for each chromosome size of maximum observed fitness scores for each permutation.}
#'  \item{max.order.pvals}{A vector of p-values for the maximum observed order statistics for each chromosome size.
#'   P-values are the proportion of permutation based maximum order statistics that exceed the observed maximum fitness score.}
#' }
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
#'                                block.ld.mat = block.ld.mat[1:10, 1:10])
#' ## run GA for observed data
#'
#' #observed data chromosome size 2
#' run.gadgets(pp.list, n.chromosomes = 5, chromosome.size = 2, results.dir = 'tmp_2',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  combined.res2 <- combine.islands('tmp_2', snp.annotations[ 1:10, ], pp.list, 2)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#'  #observed data chromosome size 3
#'  run.gadgets(pp.list, n.chromosomes = 5, chromosome.size = 3, results.dir = 'tmp_3',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  combined.res3 <- combine.islands('tmp_3', snp.annotations[ 1:10, ], pp.list, 2)
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
#'                                    block.ld.mat = block.ld.mat[1:10, 1:10])
#'
#' p2.list <- preprocess.genetic.data(perm.data.list[['permutation2']]$case,
#'                                    complement.genetic.data = perm.data.list[['permutation2']]$comp,
#'                                    block.ld.mat = block.ld.mat[1:10, 1:10])
#'
#' p3.list <- preprocess.genetic.data(perm.data.list[['permutation3']]$case,
#'                                    complement.genetic.data = perm.data.list[['permutation3']]$comp,
#'                                    block.ld.mat = block.ld.mat[1:10, 1:10])
#'
#' ##run GA for permuted data
#'
#' #permutation 1, chromosome size 2
#' run.gadgets(p1.list, n.chromosomes = 5, chromosome.size = 2, results.dir = 'p1_tmp_2',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  p1.combined.res2 <- combine.islands('p1_tmp_2', snp.annotations[ 1:10, ], p1.list, 2)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#' #permutation 1, chromosome size 3
#' run.gadgets(p1.list, n.chromosomes = 5, chromosome.size = 3, results.dir = 'p1_tmp_3',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  p1.combined.res3 <- combine.islands('p1_tmp_3', snp.annotations[ 1:10, ], p1.list, 2)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#' #permutation 2, chromosome size 2
#' run.gadgets(p2.list, n.chromosomes = 5, chromosome.size = 2, results.dir = 'p2_tmp_2',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  p2.combined.res2 <- combine.islands('p2_tmp_2', snp.annotations[ 1:10, ], p2.list, 2)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#' #permutation 2, chromosome size 3
#' run.gadgets(p2.list, n.chromosomes = 5, chromosome.size = 3, results.dir = 'p2_tmp_3',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  p2.combined.res3 <- combine.islands('p2_tmp_3', snp.annotations[ 1:10, ], p2.list, 2)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#' #permutation 3, chromosome size 2
#' run.gadgets(p3.list, n.chromosomes = 5, chromosome.size = 2, results.dir = 'p3_tmp_2',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  p3.combined.res2 <- combine.islands('p3_tmp_2', snp.annotations[ 1:10, ], p3.list, 2)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#' #permutation 3, chromosome size 3
#' run.gadgets(p3.list, n.chromosomes = 5, chromosome.size = 3, results.dir = 'p3_tmp_3',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  p3.combined.res3 <- combine.islands('p3_tmp_3', snp.annotations[ 1:10, ], p3.list, 2)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#'  ## create list of results
#'
#'  #chromosome size 2 results
#'  chrom2.list <- list(observed.data = combined.res2$fitness.score,
#'                     permutation.list = list(p1.combined.res2$fitness.score,
#'                                             p2.combined.res2$fitness.score,
#'                                             p3.combined.res2$fitness.score))
#'
#'  #chromosome size 3 results
#'  chrom3.list <- list(observed.data = combined.res3$fitness.score,
#'                     permutation.list = list(p1.combined.res3$fitness.score,
#'                                             p2.combined.res3$fitness.score,
#'                                             p3.combined.res3$fitness.score))
#'
#'  final.results <- list(chrom2.list, chrom3.list)
#'
#'  ## run global test
#'  global.test.res <- global.test(final.results, 1)
#'
#'  lapply(c('tmp_2', 'tmp_3', 'p1_tmp_2', 'p2_tmp_2', 'p3_tmp_2',
#'           'p1_tmp_3', 'p2_tmp_3', 'p3_tmp_3'), unlink, recursive = TRUE)
#'
#'
#' @export

global.test <- function(results.list, n.top.scores = 30) {

    # loop over chromosome sizes
    chrom.size.ranks <- do.call(cbind, lapply(results.list, function(chrom.size.res) {

        # error checking
        if (!"observed.data" %in% names(chrom.size.res)) {

            stop("Each element of results.list must be a list containing an element titled observed.data.")
        }

        if (!"permutation.list" %in% names(chrom.size.res)) {

            stop("Each element of results.list must be a list containing an element titled permutation.list.")
        }

        # grab the observed data

        obs <- chrom.size.res$observed.data
        if (length(obs) < n.top.scores){

            stop("n.top.scores must be <= number of unique chromosomes")

        }
        obs.score <- sum(sort(obs, decreasing = TRUE)[seq_len(n.top.scores)])

        # grab permuted data
        perm.scores <- unlist(lapply(chrom.size.res$permutation.list, function(perm.scores){

            if (length(perm.scores) < n.top.scores){

                stop("n.top.scores must be <= number of unique chromosomes")

            }
            sum(sort(perm.scores, decreasing = TRUE)[seq_len(n.top.scores)])

        }))

        # get ranks of score sums
        all.scores <- c(obs.score, perm.scores)
        score.ranks <- rank(all.scores, ties.method = "max")
        n.ranks <- length(all.scores)
        unif.ranks <- (n.ranks - score.ranks + 0.5)/n.ranks

        # return transformed scores
        return(data.frame(score = -2*log(unif.ranks)))

    }))
    colnames(chrom.size.ranks) <- NULL
    rownames(chrom.size.ranks) <- NULL

    #sum scores across chromosome sizes
    global.scores <- rowSums(chrom.size.ranks)

    # pval
    obs.test.stat <- global.scores[1]
    perm.test.stats <- global.scores[-1]
    global.pval <- sum(perm.test.stats >= obs.test.stat)/(length(perm.test.stats) + 1)

    # also look at element-wise results
    marginal.pvals <- vapply(seq_len(ncol(chrom.size.ranks)), function(x){

        chrom.size.res <- chrom.size.ranks[, x]
        obs.test.stat <- chrom.size.res[1]
        perm.test.stats <- chrom.size.res[-1]
        pval <- sum(perm.test.stats >= obs.test.stat)/(length(perm.test.stats) + 1)
        pval

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
        pval <- sum(max.perms >= max.obs)/(length(max.perms) + 1)
        pval

    }, 1.0)

    # return results list
    res.list <- list(obs.test.stat = obs.test.stat,  perm.test.stats = perm.test.stats, pval = global.pval,
                     obs.marginal.test.stats = chrom.size.ranks[1, ], perm.marginal.test.stats.mat = chrom.size.ranks[-1, ],
                     marginal.pvals = marginal.pvals,
                     max.obs.fitness = max.obs.fitness, max.perm.fitness = max.perm.fitness,
                     max.order.pvals = max.order.pvals)
    return(res.list)

}
