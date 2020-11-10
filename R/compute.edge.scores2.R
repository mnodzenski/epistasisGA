#' A function to compute edge scores for network plots of results.
#'
#' This function returns a data.table of edge weights for use in network plots of GA results.
#'
#' @param results.list A list of length d, where d is the number of chromosome sizes to be included in the network plot.
#'  Each element of the list must itself be a list whose first element \code{observed.data} is a vector of fitness scores from the
#'  the \code{unique.results} chromosome results from \code{combine.islands} for a given chromosome size. The second element \code{permutation.list}
#'  is a list containing vectors of all permutation results fitness scores, again using the \code{unique.results} results output by
#'  \code{combine.islands} for each permutation. Each data.table in the list should be subset to the top \code{n.top.scores} scores,
#'  otherwise an error will be returned.
#' @param pp.list The list output by \code{preprocess.genetic.data} run on the observed data.
#' @param n.top.scores The number of top scoring chromosomes to be used in calculating the edge.scores. Defaults to 50.
#' @param score.type A character string specifying the method for scoring edges, with options
#' 'max', 'sum', or 'logsum'. The default is 'max', but 'logsum' may also be particularly useful.
#'  Note that "logsum" is actually the log of one plus the sum of the fitness scores to avoid nodes or edges having negative
#'  weights.
#' @param epi.test.permutes The number of permutes used to compute the epistasis test p-values.
#' @param bp.param The \code{bp.param} argument to be passed to \code{run.epi.test}.
#'  If using a cluster computer, this parameter needs to be set with care. See \code{BiocParallel::bplapply} for more details
#' @return A data.table where the first two columns represent SNPs and the third column (edge.score)
#' is the edge score of a chromosome containing those SNPs.
#'
#'@examples
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
#'  chrom2.list <- list(observed.data = combined.res2$unique.results[1:3, ],
#'                     permutation.list = list(p1.combined.res2$unique.results[1:3, ],
#'                                             p2.combined.res2$unique.results[1:3, ],
#'                                             p3.combined.res2$unique.results[1:3, ]))
#'
#'  #chromosome size 3 results
#'  chrom3.list <- list(observed.data = combined.res3$unique.results[1:3, ],
#'                     permutation.list = list(p1.combined.res3$unique.results[1:3, ],
#'                                             p2.combined.res3$unique.results[1:3, ],
#'                                             p3.combined.res3$unique.results[1:3, ]))
#'
#'  final.results <- list(chrom2.list, chrom3.list)
#'
#'  ## compute edge scores
#'  edge.dt <- compute.edge.scores2(final.results, pp.list, 3)
#'
#'  lapply(c('tmp_2', 'tmp_3', 'p1_tmp_2', 'p2_tmp_2', 'p3_tmp_2',
#'           'p1_tmp_3', 'p2_tmp_3', 'p3_tmp_3'), unlink, recursive = TRUE)
#'
#' @importFrom data.table data.table rbindlist setorder
#' @importFrom matrixStats colSds
#' @importFrom utils combn
#' @export

compute.edge.scores2 <- function(results.list, pp.list, n.top.chroms = 50, score.type = "logsum",
                                 epi.test.permutes = 100, bp.param = SerialParam()) {

    ## make sure we have the correct number of chromosomes in each element of the results list
    n.chroms.vec <- lapply(results.list, function(chrom.size.list){

        n.obs.chroms <- sum(!is.na(chrom.size.list$observed.data$fitness.score))
        if (n.obs.chroms != n.top.chroms){

            stop("Each element of results.list must contain n.top.scores fitness scores.")

        }
        n.perm.chroms <- lapply(chrom.size.list$permutation.list, function(perm){

            n.perm.chroms <- sum(!is.na(perm$fitness.score))
            if (n.perm.chroms != n.top.chroms){

                stop("Each element of results.list must contain n.top.scores fitness scores.")

            }
        })

    })

    ## count the number of permutes
    n.permutes <- length(results.list[[1]]$permutation.list)

    ## compute re-scaled fitness score based on ranks compared to permutes
    rescaled.results <- lapply(results.list, function(chrom.size.list){

        obs.res <- chrom.size.list$observed.data[seq_len(n.top.chroms), ]

        # compute permutation epistasis p-values for each of 1:n.top.chroms
        # return 0.5 if we can't do the epistasis test due to all chroms being on
        # same biological chromosome
        chrom.size <- sum(grepl("snp", colnames(obs.res)))/5
        epi.pvals <- vapply(seq_len(n.top.chroms), function(x){

            chrom <- as.vector(t(obs.res[x, seq_len(chrom.size)]))
            epi.test.pval <- tryCatch(run.epi.test(chrom, pp.list, n.permutes = epi.test.permutes, bp.param = bp.param)$pval,
                                     error = function(e) 0.5)
            return(epi.test.pval)

        }, 1.0)

        obs.scores <- matrix(rep(sort(obs.res$fitness.score, decreasing = TRUE), n.permutes),
                             nrow = n.permutes, byrow = TRUE)
        perm.scores <- do.call(rbind, lapply(chrom.size.list$permutation.list, function(x) sort(x$fitness.score[seq_len(n.top.chroms)], decreasing = TRUE)))
        rescaled.scores <- colSums(perm.scores >= obs.scores)/(n.permutes + 1)

        ## set zeros to small fraction
        if (any(rescaled.scores == 0)){
            rescaled.scores[rescaled.scores == 0] <- 1/(n.permutes + 1)
        }

        ## set 1's to large fraction
        if (any(rescaled.scores == 1)){

            rescaled.scores[rescaled.scores == 1] <- (n.permutes)/(n.permutes + 1)

        }

        ## take -2 log
        rescaled.scores <- -2*log(rescaled.scores)

        ## updated the observed results data.table
        obs.res$fitness.score <- rescaled.scores
        obs.res$chrom.weight <- 1 - epi.pvals
        return(obs.res)

    })

    all.edge.weights <- rbindlist(lapply(seq_along(rescaled.results), function(d){

        chrom.size.res <- rescaled.results[[d]]
        chrom.size <- sum(grepl("snp", colnames(chrom.size.res)))/5
        rbindlist(lapply(seq_len(n.top.chroms), function(res.row) {

            chrom.res <- chrom.size.res[res.row, ]
            fs <- chrom.res$fitness.score
            chrom.weight <- chrom.res$chrom.weight
            these.cols <- seq_len(chrom.size)
            chrom <- as.vector(t(chrom.res[, ..these.cols]))
            chrom.pairs <- data.table(t(combn(chrom, 2)))
            chrom.pairs[ , `:=`(raw.fitness.score = fs, chrom.weight = chrom.weight,
                                weighted.fs = chrom.weight*fs)]

            #take weighted average score for each pair
            out.dt <- chrom.pairs[ , list(fitness.score = sum(weighted.fs)/sum(chrom.weight)),
                                   list(V1, V2)]
            return(out.dt)

        }))

    }))

    # compute edge score based on score.type
    if (score.type == "max"){

        out.dt <- all.edge.weights[ , list(edge.score = max(fitness.score)), list(V1, V2)]
        setorder(out.dt, -edge.score)

    } else if (score.type == "sum"){

        out.dt <- all.edge.weights[ , list(edge.score = sum(fitness.score)/length(results.list)),
                                    list(V1, V2)]
        setorder(out.dt, -edge.score)

    } else if (score.type == "logsum"){

        out.dt <- all.edge.weights[ , list(edge.score = log(1 + (sum(fitness.score)/length(results.list)))),
                                    list(V1, V2)]
        setorder(out.dt, -edge.score)

    }

    return(out.dt)

}


