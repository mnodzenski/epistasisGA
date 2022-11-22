#' A function to run a global test of the null hypothesis that there are no joint SNP-disease associations
#' across a range of chromosome sizes
#'
#' This function runs a global test of the null hypothesis that there are no joint SNP-disease associations
#' across a range of chromosome sizes
#'
#' @param results.list A list of length d, where d is the number of chromosome
#' sizes to be included in a global test. Each element of the list must itself
#' be a list whose first element \code{observed.data} is a vector of fitness
#' scores from \code{combine.islands} for the observed data for a given
#' chromosome size. The second element \code{permutation.list} is a list
#' containing vectors of all permutation results fitness scores, again using the
#' results output by \code{combine.islands} for each permutation.
#' @param n.top.scores The number of top scoring chromosomes, for each
#' chromosome size, to be used in calculating the global test. Defaults to 10.
#' @return A list containing the following:
#' \describe{
#'  \item{obs.test.stat}{The observed test statistic.}
#'  \item{perm.test.stats}{A vector of test statistics from permuted data.}
#'  \item{pval}{The p-value for the global test.}
#'  \item{obs.marginal.test.stats}{A vector of observed test statistics for
#'  each chromosome size.}
#'  \item{perm.marginal.test.stats.mat}{A matrix of test statistics for the
#'  permutation datasets, where rows correspond to
#'  permutations and columns correspond to chromosome sizes.}
#'  \item{marginal.pvals}{A vector containing marignal p-values for each
#'  chromosome size.}
#'  \item{max.obs.fitness}{A vector of the maximum fitness score for each
#'  chromosome size in the observed data.}
#'  \item{max.perm.fitness}{A list of vectors for each chromosome size of
#'  maximum observed fitness scores for each permutation.}
#'  \item{max.order.pvals}{A vector of p-values for the maximum observed order
#'  statistics for each chromosome size. P-values are the proportion of
#'  permutation based maximum order statistics that exceed the observed maximum
#'  fitness score.}
#'  \item{boxplot.grob}{A grob of a ggplot plot of the observed vs permuted
#'  fitness score densities for each chromosome size.}
#'  \item{chrom.size.k}{A vector indicating the number of top scores (k)
#'  from each chromosome size that the test used.
#'  This will be equal to \code{n.top.scores} unless GADGETS returns fewer than
#'  \code{n.top.scores} unique chromosomes for
#'  the observed data or any permute, in which case the chromosome size-specific
#'   value will be equal to the smallest number of unique chromosomes returned.}
#'  \item{max.perm.95th.pctl}{The 95th percentile of the permutation maximum
#'  order statistics for each chromosome size.}
#' }
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' case <- as.matrix(case)
#' dad <- as.matrix(dad)
#' mom <- as.matrix(mom)
#' data(snp.annotations)
#' set.seed(1400)
#'
# #preprocess data
#' pp.list <- preprocess.genetic.data(case[, 1:10],
#'                                father.genetic.data = dad[ , 1:10],
#'                                mother.genetic.data = mom[ , 1:10],
#'                                ld.block.vec = c(10))
#' ## run GA for observed data
#'
#' #observed data chromosome size 2
#' run.gadgets(pp.list, n.chromosomes = 5, chromosome.size = 2,
#'        results.dir = 'tmp_2',
#'        cluster.type = 'interactive',
#'        registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  combined.res2 <- combine.islands('tmp_2', snp.annotations[ 1:10, ],
#'                                   pp.list, 2)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#'  #observed data chromosome size 3
#'  run.gadgets(pp.list, n.chromosomes = 5, chromosome.size = 3,
#'        results.dir = 'tmp_3',
#'        cluster.type = 'interactive',
#'        registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  combined.res3 <- combine.islands('tmp_3', snp.annotations[ 1:10, ],
#'                                   pp.list, 2)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#' #create three permuted datasets
#' set.seed(1400)
#' perm.data.list <- permute.dataset(pp.list, "perm_data",
#'                                   n.permutations = 3)
#'
#' #pre-process permuted data
#' case.p1 <- readRDS("perm_data/case.permute1.rds")
#' comp.p1 <- readRDS("perm_data/complement.permute1.rds")
#' p1.list <- preprocess.genetic.data(case.p1,
#'                                    complement.genetic.data = comp.p1,
#'                                     ld.block.vec = c(10))
#'
#' case.p2 <- readRDS("perm_data/case.permute2.rds")
#' comp.p2 <- readRDS("perm_data/complement.permute2.rds")
#' p2.list <- preprocess.genetic.data(case.p2,
#'                                    complement.genetic.data = comp.p2,
#'                                     ld.block.vec = c(10))
#'
#' case.p3 <- readRDS("perm_data/case.permute3.rds")
#' comp.p3 <- readRDS("perm_data/complement.permute3.rds")
#' p3.list <- preprocess.genetic.data(case.p3,
#'                                    complement.genetic.data = comp.p3,
#'                                    ld.block.vec = c(10))
#'
#' ##run GA for permuted data
#'
#' #permutation 1, chromosome size 2
#' run.gadgets(p1.list, n.chromosomes = 5, chromosome.size = 2,
#'        results.dir = 'p1_tmp_2',
#'        cluster.type = 'interactive',
#'        registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  p1.combined.res2 <- combine.islands('p1_tmp_2', snp.annotations[ 1:10, ],
#'                                      p1.list, 2)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#' #permutation 1, chromosome size 3
#' run.gadgets(p1.list, n.chromosomes = 5, chromosome.size = 3,
#'        results.dir = 'p1_tmp_3',
#'        cluster.type = 'interactive',
#'        registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  p1.combined.res3 <- combine.islands('p1_tmp_3', snp.annotations[ 1:10, ],
#'                                      p1.list, 2)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#' #permutation 2, chromosome size 2
#' run.gadgets(p2.list, n.chromosomes = 5, chromosome.size = 2,
#'        results.dir = 'p2_tmp_2',
#'        cluster.type = 'interactive',
#'        registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  p2.combined.res2 <- combine.islands('p2_tmp_2', snp.annotations[ 1:10, ],
#'                                      p2.list, 2)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#' #permutation 2, chromosome size 3
#' run.gadgets(p2.list, n.chromosomes = 5, chromosome.size = 3,
#'        results.dir = 'p2_tmp_3',
#'        cluster.type = 'interactive',
#'        registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  p2.combined.res3 <- combine.islands('p2_tmp_3', snp.annotations[ 1:10, ],
#'                                      p2.list, 2)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#' #permutation 3, chromosome size 2
#' run.gadgets(p3.list, n.chromosomes = 5, chromosome.size = 2,
#'        results.dir = 'p3_tmp_2',
#'        cluster.type = 'interactive',
#'        registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  p3.combined.res2 <- combine.islands('p3_tmp_2', snp.annotations[ 1:10, ],
#'                                      p3.list, 2)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#' #permutation 3, chromosome size 3
#' run.gadgets(p3.list, n.chromosomes = 5, chromosome.size = 3,
#'        results.dir = 'p3_tmp_3',
#'        cluster.type = 'interactive',
#'        registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  p3.combined.res3 <- combine.islands('p3_tmp_3', snp.annotations[ 1:10, ],
#'                                      p3.list, 2)
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
#'           'p1_tmp_3', 'p2_tmp_3', 'p3_tmp_3', 'perm_data'), unlink,
#'           recursive = TRUE)
#'
#'
#' @importFrom data.table data.table rbindlist
#' @importFrom ggplot2 ggplot aes facet_wrap geom_boxplot ggplotGrob
#' @importFrom grid grid.draw
#' @importFrom stats quantile
#' @export

global.test <- function(results.list, n.top.scores = 10) {

    # decide on the n.top.scores per chromosome size
    chrom.size.n.top.scores <- rep(n.top.scores, length(results.list))
    for (i in seq_along(chrom.size.n.top.scores)){

        chrom.size.res <- results.list[[i]]
        # error checking
        if (!"observed.data" %in% names(chrom.size.res)) {

            stop("Each element of results.list must be a list containing an element titled observed.data.")
        }

        if (!"permutation.list" %in% names(chrom.size.res)) {

            stop("Each element of results.list must be a list containing an element titled permutation.list.")
        }

        # re-set n.top.chroms if needed
        n.obs.scores <- length(chrom.size.res$observed.data)
        n.perm.scores <- vapply(chrom.size.res$permutation.list, length, 1)
        min.n.scores <- min(c(n.obs.scores, n.perm.scores))
        if (min.n.scores < n.top.scores){

            chrom.size.n.top.scores[i] <- min.n.scores
            if (n.obs.scores == min.n.scores){

                message(paste("only", min.n.scores, "observed.data fitness scores for element",
                              i, "of results list, re-setting n.top.chroms to",
                              min.n.scores))

            } else if (any(n.perm.scores == min.n.scores)){

                message(paste("only", min.n.scores, "permutation.list fitness scores for element",
                              i, "of results list, re-setting n.top.chroms to",
                              min.n.scores))

            }

        }

    }

    # loop over chromosome sizes
    chrom.size.ranks <- do.call(cbind, lapply(seq_along(results.list),
                                              function(x) {

        # get the chromosome size specific data
        chrom.size.res <- results.list[[x]]
        n.top <- chrom.size.n.top.scores[x]

        # grab the observed data
        obs <- chrom.size.res$observed.data
        obs.score <- sum(sort(obs, decreasing = TRUE)[seq_len(n.top)])

        # grab permuted data
        perm.scores <- unlist(lapply(chrom.size.res$permutation.list,
                                     function(perm.scores){

            sum(sort(perm.scores, decreasing = TRUE)[seq_len(n.top)])

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
    B <- sum(perm.test.stats >= obs.test.stat)
    N <- length(perm.test.stats) + 1
    global.pval <- (B + 1)/N

    # also look at element-wise results
    marginal.pvals <- vapply(seq_len(ncol(chrom.size.ranks)), function(x){

        chrom.size.res <- chrom.size.ranks[, x]
        obs.test.stat <- chrom.size.res[1]
        perm.test.stats <- chrom.size.res[-1]
        B <- sum(perm.test.stats >= obs.test.stat)
        N <- length(perm.test.stats) + 1
        pval <- (B + 1)/N
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

    # also return 95th percentile of the null maxima
    max.perm.95th.pctl <- apply(max.perm.fitness, 2, function(x){

        quantile(x, 0.95)

    })
    names(max.perm.95th.pctl) <- NULL

    # pvals for max order statistics
    max.order.pvals <- vapply(seq_along(max.obs.fitness), function(chrom.size) {

        max.obs <- max.obs.fitness[chrom.size]
        max.perms <- max.perm.fitness[ , chrom.size]
        B <- sum(max.perms >= max.obs)
        N <- length(max.perms) + 1
        pval <- (B + 1)/N
        pval

    }, 1.0)

    # plots of the observed vs permute distributions
    plot.data <- rbindlist(lapply(seq_along(results.list), function(x) {

        # get the chromosome size specific data
        chrom.size.res <- results.list[[x]]

        # grab the observed data
        obs <- chrom.size.res$observed.data
        obs.dt <- data.table(chrom.res.element = x, data.type = "observed",
                             fitness.score = obs)

        # grab permuted data
        perm.dt <- rbindlist(lapply(chrom.size.res$permutation.list,
                                    function(perm.scores){

            data.table(chrom.res.element = x, data.type = "permuted",
                       fitness.score = perm.scores)

        }))

        # put together data.table for plotting
        plot.dt <- rbind(obs.dt, perm.dt)
        return(plot.dt)

    }))

    boxplot.plot <- ggplot(plot.data, aes(x = data.type, y = fitness.score)) +
        geom_boxplot() + facet_wrap(chrom.res.element ~ ., scales = "free")

    # don't need copies of the global environment vars
    boxplot.plot$plot_env <- new.env()

    #store as grob
    boxplot.grob <- ggplotGrob(boxplot.plot)

    # return results list
    res.list <- list(obs.test.stat = obs.test.stat,
                     perm.test.stats = perm.test.stats, pval = global.pval,
                     obs.marginal.test.stats = chrom.size.ranks[1, ],
                     perm.marginal.test.stats.mat = chrom.size.ranks[-1, ],
                     marginal.pvals = marginal.pvals,
                     max.obs.fitness = max.obs.fitness,
                     max.perm.fitness = max.perm.fitness,
                     max.order.pvals = max.order.pvals,
                     boxplot.grob = boxplot.grob,
                     chrom.size.k = chrom.size.n.top.scores,
                     max.perm.95th.pctl = max.perm.95th.pctl)
    return(res.list)

}
