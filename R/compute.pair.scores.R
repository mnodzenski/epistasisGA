#' A function to compute SNP-pair scores for network plots of results.
#'
#' This function returns a data.table of graphical SNP-pair scores for use in network plots of GADGETS results.
#'
#' @param results.list A list of length d, where d is the number of chromosome sizes to be included in the network plot.
#'  Each element of the list must be a data.table from \code{combine.islands} for a given chromosome size.
#'  Each data.table in the list should be subset to the top \code{n.top.scores} scores,
#'  otherwise an error will be returned.
#' @param pp.list The list output by \code{preprocess.genetic.data} run on the observed data.
#' @param n.top.chroms The number of top scoring chromosomes to be used in calculating the edge.scores. Defaults to 50.
#' @param score.type A character string specifying the method for aggregating SNP-pair scores across chromosome sizes. Options are
#' 'max', 'sum', or 'logsum', defaulting to "logsum". For a given SNP-pair, it's graphical score will be the \code{score.type} of all
#' graphical scores of chromosomes containing that pair across chromosome sizes. The choice of 'logsum' rather than 'sum'
#' may be useful in cases where there are multiple risk-sets, and one is found much more frequently. However, it may be of interest to examine
#' plots using both \code{score.type} approaches. Note that "logsum" is actually the log of one plus the sum of the SNP-pair scores to avoid nodes or
#' edges having negative weights.
#' @param epi.test.permutes The number of permutes used to compute the epistasis test p-values. Defaults to 100.
#' @param bp.param The \code{bp.param} argument to be passed to \code{run.epi.test}.
#'  If using a cluster computer, this parameter needs to be set with care. See \code{BiocParallel::bplapply} for more details
#' @param pval.thresh A numeric value between 0 and 1 specifying the epistasis test p-value threshold for a
#' chromosome to contribute to the network. Any chromosomes with epistasis p-value greater than \code{pval.thresh}
#' will not contribute to network plots. The argument defaults to 0.05.
#' @return A data.table where the first four columns represent SNPs and the fifth column (edge.score)
#' is the graphical SNP-pair score.
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
#' #preprocess data
#' target.snps <- c(1:3, 30:32, 60:62, 85)
#' pp.list <- preprocess.genetic.data(case[, target.snps], father.genetic.data = dad[ , target.snps],
#'                                mother.genetic.data = mom[ , target.snps],
#'                                block.ld.mat = block.ld.mat[target.snps , target.snps])
#' ## run GA for observed data
#'
#' #observed data chromosome size 2
#' run.gadgets(pp.list, n.chromosomes = 5, chromosome.size = 2, results.dir = 'tmp_2',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3,
#'        n.migrations = 0)
#'  combined.res2 <- combine.islands('tmp_2', snp.annotations[ target.snps, ], pp.list)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#'  #observed data chromosome size 3
#'  run.gadgets(pp.list, n.chromosomes = 5, chromosome.size = 3, results.dir = 'tmp_3',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3,
#'        n.migrations = 0)
#'  combined.res3 <- combine.islands('tmp_3', snp.annotations[ target.snps, ], pp.list)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#'  ## create list of results
#'
#'  final.results <- list(combined.res2[1:3, ], combined.res3[1:3, ])
#'
#'  ## compute edge scores
#'  edge.dt <- compute.pair.scores(final.results, pp.list, 3, pval.thresh = 1)
#'
#'  lapply(c('tmp_2', 'tmp_3'), unlink, recursive = TRUE)
#'
#' @importFrom data.table data.table rbindlist setorder
#' @importFrom matrixStats colSds
#' @importFrom utils combn
#' @importFrom BiocParallel SerialParam
#' @export

compute.pair.scores <- function(results.list, pp.list, n.top.chroms = 50, score.type = "logsum",
                                 epi.test.permutes = 100, bp.param = SerialParam(), pval.thresh = 0.05) {

    ## make sure we have the correct number of chromosomes in each element of the results list
    n.chroms.vec <- lapply(results.list, function(chrom.size.data){

        n.obs.chroms <- sum(!is.na(chrom.size.data$fitness.score))
        if (n.obs.chroms != n.top.chroms){

            stop("Each element of results.list must contain n.top.scores fitness scores.")

        }

    })

    ## compute re-scaled fitness scores based on epistasis test
    rescaled.results <- lapply(results.list, function(obs.res){

        # compute permutation epistasis p-values for each of 1:n.top.chroms
        # return 0.5 if we can't do the epistasis test due to all chroms being on
        # same biological chromosome
        chrom.size <- sum(grepl("snp", colnames(obs.res)))/5
        epi.pvals <- vapply(seq_len(n.top.chroms), function(x){

            these.cols <- seq_len(chrom.size)
            chrom <- as.vector(t(obs.res[x, ..these.cols]))
            epi.test.pval <- tryCatch(epistasis.test(chrom, pp.list, n.permutes = epi.test.permutes, bp.param = bp.param)$pval,
                                     error = function(e) {

                                         error.message = conditionMessage(e)
                                         if (error.message == "cannot run test, all SNPs are in LD"){

                                             0.5

                                         } else {

                                             stop(e)

                                         }

                                     })
            if (isTRUE(all.equal(0, epi.test.pval))){

                epi.test.pval <- 1/(epi.test.permutes + 1)

            }
            if (isTRUE(all.equal(1, epi.test.pval))){

                epi.test.pval <- 1 - (1/(epi.test.permutes + 1))

            }
            return(epi.test.pval)

        }, 1.0)

        ## take -2 log
        rescaled.scores <- -2*log(epi.pvals)

        ## update the observed results data.table
        obs.res$fitness.score <- rescaled.scores
        obs.res <- obs.res[epi.pvals <= pval.thresh, ]
        return(obs.res)

    })

    #make sure we have some edges
    n.edges <- vapply(rescaled.results, nrow, 1)
    if (sum(n.edges) == 0){
        stop("No SNP pairs meet p-value threshold")
    }

    all.edge.weights <- rbindlist(lapply(seq_along(rescaled.results), function(d){

        chrom.size.res <- rescaled.results[[d]]
        chrom.size <- sum(grepl("snp", colnames(chrom.size.res)))/5
        rbindlist(lapply(seq_len(n.top.chroms), function(res.row) {

            chrom.res <- chrom.size.res[res.row, ]
            fs <- chrom.res$fitness.score
            these.cols <- seq_len(chrom.size)
            chrom <- as.vector(t(chrom.res[, ..these.cols]))
            chrom.pairs <- data.table(t(combn(chrom, 2)))
            colnames(chrom.pairs) <- c("SNP1", "SNP2")
            rsid.cols <- seq(chrom.size + 1, 2*chrom.size)
            rsids <- as.vector(t(chrom.res[, ..rsid.cols]))
            if (chrom.size == 2){

                rsid.dt <- data.table(t(apply(chrom.pairs, 2, function(x) rsids[match(x, chrom)])))

            } else {

                rsid.dt <- data.table(apply(chrom.pairs, 2, function(x) rsids[match(x, chrom)]))

            }
            colnames(rsid.dt) <- c("SNP1.rsid", "SNP2.rsid")
            chrom.pairs <- cbind(chrom.pairs, rsid.dt)
            chrom.pairs[ , `:=`(raw.fitness.score = fs)]

            #take average score for each pair
            out.dt <- chrom.pairs[ , list(fitness.score = mean(raw.fitness.score)),
                                   list(SNP1, SNP2, SNP1.rsid, SNP2.rsid)]
            return(out.dt)

        }))

    }))

    # remove the NA's
    all.edge.weights <- all.edge.weights[!is.na(all.edge.weights$fitness.score), ]

    # compute edge score based on score.type
    if (score.type == "max"){

        out.dt <- all.edge.weights[ , list(edge.score = max(fitness.score)), list(SNP1, SNP2, SNP1.rsid, SNP2.rsid)]
        setorder(out.dt, -edge.score)

    } else if (score.type == "sum"){

        out.dt <- all.edge.weights[ , list(edge.score = sum(fitness.score)/length(results.list)),
                                    list(SNP1, SNP2, SNP1.rsid, SNP2.rsid)]
        setorder(out.dt, -edge.score)

    } else if (score.type == "logsum"){

        out.dt <- all.edge.weights[ , list(edge.score = log(1 + (sum(fitness.score)/length(results.list)))),
                                    list(SNP1, SNP2, SNP1.rsid, SNP2.rsid)]
        setorder(out.dt, -edge.score)

    }

    return(out.dt)

}


