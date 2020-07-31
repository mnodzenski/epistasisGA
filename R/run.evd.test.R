#' A function to run tests of the null hypothesis that observed fitness scores
#' were not drawn from the null extreme value distribution of maximum fitness scores,
#' using permuted null datasets to estimate the EVD parameters.
#'
#' This function runs tests of the null hypothesis that observed fitness scores
#' were not drawn from the null extreme value distribution of maximum fitness scores,
#' using permuted null datasets to estimate the EVD parameters.
#'
#' @param results.list A list of length d, where d is the number of chromosome sizes to be included.
#'  Each element of the list must itself be a list whose first element \code{observed.data} is a vector containing
#'  the fitness scores from the \code{unique.results} chromosome results from \code{combine.islands} for a given chromosome size.
#'  The second element \code{permutation.list} is a list containing all permutation results fitness scores, again using the
#'  \code{unique.results} results output by \code{combine.islands} for each permutation.
#' @param ... Additional arguments to be passed to \code{evd::fgev}.
#' @return A list of length d, with each element corresponding to a chromosome size, and
#' each sub-element containing a list with the following:
#' \describe{
#'  \item{pval}{A vector of the upper tail probabilities for each of the input observed fitness scores
#'  based on the null extreme value distribution estimated using the permuted data.}
#'  \item{null.evd.fit}{The fitted evd object returned by \code{fgev} run on the maxima of the permuted
#'  data fitness scores.}
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
#'  chrom2.list <- list(observed.data = combined.res2$unique.results$fitness.score,
#'                     permutation.list = list(p1.combined.res2$unique.results$fitness.score,
#'                                             p2.combined.res2$unique.results$fitness.score,
#'                                             p3.combined.res2$unique.results$fitness.score))
#'
#'  #chromosome size 3 results
#'  chrom3.list <- list(observed.data = combined.res3$unique.results$fitness.score,
#'                     permutation.list = list(p1.combined.res3$unique.results$fitness.score,
#'                                             p2.combined.res3$unique.results$fitness.score,
#'                                             p3.combined.res3$unique.results$fitness.score))
#'
#'  final.results <- list(chrom2.list, chrom3.list)
#'
#'  ## run test
#'  evd.test.res <- run.evd.test(final.results)
#'
#'  lapply(c('tmp_2', 'tmp_3', 'p1_tmp_2', 'p2_tmp_2', 'p3_tmp_2',
#'           'p1_tmp_3', 'p2_tmp_3', 'p3_tmp_3'), unlink, recursive = TRUE)
#'
#'
#' @importFrom matrixStats rowMaxs
#' @importFrom stats ecdf
#' @importFrom evd fgev pgev
#' @export

run.evd.test <- function(results.list, ...) {

    # loop over chromosome sizes
    chrom.size.list <- lapply(results.list, function(chrom.size.res) {

        # error checking
        if (!"observed.data" %in% names(chrom.size.res)) {

            stop("Each element of results.list must be a list containing an element titled observed.data.")
        }

        if (!"permutation.list" %in% names(chrom.size.res)) {

            stop("Each element of results.list must be a list containing an element titled permutation.list.")
        }

        # grab the observed data
        obs.fitness.scores <- chrom.size.res$observed.data

        # grab the max fitness score from each permuted dataset
        perm.list <- chrom.size.res$permutation.list
        perm.maxes <- vapply(perm.list, max, 1.0)

        # estimate the null evd
        null.evd <- fgev(perm.maxes, ...)
        null.loc <- null.evd$param['loc']
        null.scale <- null.evd$param['scale']
        null.shape <- null.evd$param['shape']

        # compute upper tail probs for observed fitness scores
        pvals <- pgev(obs.fitness.scores, loc = null.loc,
                      scale = null.scale, shape = null.shape,
                      lower.tail = FALSE)

        # return pvals, fitted evd object
        list(pval = pvals, null.evd.fit = null.evd)

    })

    #return final list
    return(chrom.size.list)

}
