#' A function to compute SNP-pair scores for network plots of results.
#'
#' This function returns a data.table of graphical SNP-pair scores for use in
#'  network plots of GADGETS results.
#'
#' @param results.list A list of length d, where d is the number of
#' chromosome sizes to be included in the network plot. Each element of the list
#' should be a data.table from \code{combine.islands} for a given chromosome
#' size. Each data.table in the list should be subset to only include those
#' chromosomes whose fitness scores are high enough to contribute to the network
#' plot. The selection of the chromosomes that contribute to these plots
#' is at the analyst's discretion. We have found success in just using the
#' top 10 scoring chromosomes, and also by restricting attention to those
#' chromosomes that exceed the 95th percentile of the maxima observed after
#' running GADGETS on data-sets permuted under a global, no-association null.
#' For the latter, see also function \code{global.test}.
#' @param preprocessed.list The list output by \code{preprocess.genetic.data}
#' run on the observed data.
#' @param score.type A character string specifying the method for aggregating
#' SNP-pair scores across chromosome sizes. Options are 'max', 'sum', or
#' 'logsum', defaulting to 'logsum'. For a given SNP-pair, it's graphical score
#' will be the \code{score.type} of all graphical scores of chromosomes
#' containing that pair across chromosome sizes. The choice of 'logsum' rather
#' than 'sum' may be useful in cases where there are multiple risk-sets, and one
#' is found much more frequently. However, it may be of interest to examine
#' plots using both \code{score.type} approaches.
#' @param pval.thresh A numeric value between 0 and 1 specifying the epistasis
#' test p-value threshold for a chromosome to contribute to the network. Any
#' chromosomes with epistasis p-value greater than \code{pval.thresh}
#' will not contribute to network plots. The argument defaults to 0.05.
#' It must be <= 0.6.
#' @param n.permutes The number of permutations on which to base the epistasis
#' tests. Defaults to 10000.
#' @param n.different.snps.weight The number by which the number of different
#' SNPs between a case and complement/unaffected sibling
#' is multiplied in computing the family weights. Defaults to 2.
#' @param n.both.one.weight The number by which the number of SNPs equal to 1
#' in both the case and complement/unaffected sibling
#' is multiplied in computing the family weights. Defaults to 1.
#' @param weight.function.int An integer used to assign family weights.
#' Specifically, we use \code{weight.function.int} in a function that takes
#' the weighted sum of the number of different SNPs and SNPs both equal to one
#' as an argument, denoted as x, and
#' returns a family weight equal to \code{weight.function.int}^x. Defaults to 2.
#' @param recessive.ref.prop The proportion to which the observed proportion of
#' informative cases with the provisional risk genotype(s) will be compared
#' to determine whether to recode the SNP as recessive. Defaults to 0.75.
#' @param recode.test.stat For a given SNP, the minimum test statistic required
#' to recode and recompute the fitness score using recessive coding.
#' Defaults to 1.64.
#' @param bp.param The BPPARAM argument to be passed to bplapply.
#' See \code{BiocParallel::bplapply} for more details.
#' @param null.mean.vec.list (experimental) A list, equal in length to
#' \code{results.list}, where the i^th element of the list is the vector of null
#' means (stored in the 'null.mean.sd.info.rds') corresponding to the d
#' (chromosome size) used to generate the results stored in the i^th
#' element of \code{results.list}. This only needs to be specified if based on
#' the experimental E-GADGETS method, and otherwise can be left at its default.
#' @param null.sd.vec.list (experimental) A list, equal in length to
#' \code{results.list}, where the i^th element of the list is the vector of null
#' standard deviations (stored in the 'null.mean.sd.info.rds') corresponding to
#' the d (chromosome size) used to generate the results stored in the i^th
#' element of \code{results.list}. This only needs to be specified if based on
#' the experimental E-GADGETS method, and otherwise can be left at its default.
#' @return A list of two elements:
#' \describe{
#'  \item{pair.scores}{A data.table containing SNP-pair graphical scores,
#'  where the first four columns represent SNPs and the fifth column (pair.score)
#' is the graphical SNP-pair score.}
#'  \item{snp.scores}{A data.table containing individual SNP graphical scores,
#'  where the first two columns represent SNPs and the third column (snp.score)
#'  is the graphical SNP score.}
#' }
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' data(snp.annotations)
#' set.seed(1400)
#'
#' # preprocess data
#' target.snps <- c(1:3, 30:32, 60:62, 85)
#' preprocessed.list <- preprocess.genetic.data(as.matrix(case[, target.snps]),
#'                         father.genetic.data = as.matrix(dad[ , target.snps]),
#'                         mother.genetic.data = as.matrix(mom[ , target.snps]),
#'                      ld.block.vec = c(3, 3, 3, 1))
#' ## run GA for observed data
#'
#' #observed data chromosome size 2
#' run.gadgets(preprocessed.list, n.chromosomes = 5, chromosome.size = 2,
#'        results.dir = 'tmp_2',
#'        cluster.type = 'interactive',
#'        registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  combined.res2 <- combine.islands('tmp_2',
#'                      snp.annotations[ target.snps, ], preprocessed.list, 2)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#'  #observed data chromosome size 3
#'  run.gadgets(preprocessed.list, n.chromosomes = 5,
#'        chromosome.size = 3, results.dir = 'tmp_3',
#'        cluster.type = 'interactive',
#'        registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  combined.res3 <- combine.islands('tmp_3', snp.annotations[ target.snps, ],
#'                                   preprocessed.list, 2)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#' ## create list of results
#'
#' final.results <- list(combined.res2[1:3, ], combined.res3[1:3, ])
#'
#'  ## compute edge scores
#'  edge.dt <- compute.graphical.scores(final.results,
#'                                      preprocessed.list,
#'                                      pval.thresh = 0.5)
#'
#' lapply(c("tmp_2", "tmp_3"), unlink, recursive = TRUE)
#'
#' @importFrom data.table data.table rbindlist setorder
#' @importFrom matrixStats colSds
#' @importFrom utils combn
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom data.table melt
#' @export compute.graphical.scores

compute.graphical.scores <- function(results.list, preprocessed.list,
                                     score.type = "logsum", pval.thresh = 0.05,
                                     n.permutes = 10000,
                                     n.different.snps.weight = 2,
                                     n.both.one.weight = 1,
                                     weight.function.int = 2,
                                     recessive.ref.prop = 0.75,
                                     recode.test.stat = 1.64,
                                     bp.param = bpparam(),
                                     null.mean.vec.list = NULL,
                                     null.sd.vec.list = NULL) {

    if (pval.thresh > 0.6){

        stop("pval.thresh must be <= 0.6")
    }

    #need to specify both or neither of null mean and sd vec
    if (is.null(null.mean.vec.list) & !is.null(null.sd.vec.list) |
        !is.null(null.mean.vec.list) & is.null(null.sd.vec.list)){

        stop("null.mean.vec.list and null.sd.vec.list must both be NULL or both not NULL")

    }

    #make sure mean and sd vec lists are the same length as results list
    if (!is.null(null.mean.vec.list)){

        if (length(null.mean.vec.list) != length(results.list) |
            length(null.sd.vec.list) != length(results.list)){

            stop("null.mean.vec.list and null.sd.vec.list must have the same length as results list.")

        }

    }

    ## divide up results into list of chromosome lists
    chrom.list <- lapply(results.list, function(chrom.size.data){

        n.obs.chroms <- sum(!is.na(chrom.size.data$fitness.score))
        chrom.size <- which(colnames(chrom.size.data) == "snp1.rsid") - 1
        these.cols <- seq_len(chrom.size)
        chrom.mat <- as.matrix(chrom.size.data[, ..these.cols])
        chrom.list <- split(chrom.mat, seq_len(nrow(chrom.mat)))
        return(chrom.list)
    })

    ## compute graphical scores based on epistasis/GxE test
    n2log.epi.pvals <- bplapply(seq_along(chrom.list),
                                function(i, chrom.list,
                                         preprocessed.list,
                                         null.mean.vec.list,
                                         null.sd.vec.list,
                                         n.permutes,
                                         n.different.snps.weight,
                                         n.both.one.weight,
                                         weight.function.int,
                                         recessive.ref.prop,
                                         recode.test.stat){

        chrom.size.list <- chrom.list[[i]]
        if (!is.null(null.mean.vec.list)){

            null.mean.vec <- null.mean.vec.list[[i]]
            null.sd.vec <- null.sd.vec.list[[i]]

        } else {

            null.mean.vec <- NULL
            null.sd.vec <- NULL

        }
        n2log_epistasis_pvals(chrom.size.list, preprocessed.list, n.permutes,
                              n.different.snps.weight, n.both.one.weight,
                              weight.function.int, recessive.ref.prop,
                              recode.test.stat, null.mean.vec, null.sd.vec)

    }, chrom.list = chrom.list, preprocessed.list = preprocessed.list,
    null.mean.vec.list = null.mean.vec.list,
    null.sd.vec.list = null.sd.vec.list, n.permutes = n.permutes,
    n.different.snps.weight = n.different.snps.weight,
    n.both.one.weight = n.both.one.weight,
    weight.function.int = weight.function.int,
    recessive.ref.prop = recessive.ref.prop,
    recode.test.stat = recode.test.stat, BPPARAM = bp.param)

   ## add those scores to the obs data
   for (i in seq_along(n2log.epi.pvals)){

       results.list[[i]]$graphical.score <- n2log.epi.pvals[[i]]
       results.list[[i]] <- results.list[[i]][
           results.list[[i]]$graphical.score >= -2*log(pval.thresh), ]

   }

    # make sure we have some edges
    n.edges <- vapply(results.list, nrow, 1)
    if (sum(n.edges) == 0) {
        stop("No SNP pairs meet p-value threshold")
    }

    # get rid of d where we have no edges
    results.list <- results.list[n.edges > 0]

    # edge weights
    edge.weights <- rbindlist(bplapply(seq_along(results.list),
                                       function(d, results.list){

        chrom.size.res <- results.list[[d]]
        n.top.chroms <- nrow(chrom.size.res)
        chrom.size <- which(colnames(chrom.size.res) == "snp1.rsid") - 1
        chrom.size.pairs <- rbindlist(lapply(seq_len(n.top.chroms),
                                             function(res.row) {

            chrom.res <- chrom.size.res[res.row, ]
            score <- chrom.res$graphical.score
            these.cols <- seq_len(chrom.size)
            chrom <- as.vector(t(chrom.res[, ..these.cols]))
            chrom.pairs <- data.table(t(combn(chrom, 2)))
            colnames(chrom.pairs) <- c("SNP1", "SNP2")
            rsid.cols <- seq(chrom.size + 1, 2*chrom.size)
            rsids <- as.vector(t(chrom.res[, ..rsid.cols]))
            if (chrom.size == 2){

                rsid.dt <- data.table(t(apply(chrom.pairs, 2, function(x)
                    rsids[match(x, chrom)])))

            } else {

                rsid.dt <- data.table(apply(chrom.pairs, 2, function(x)
                    rsids[match(x, chrom)]))

            }
            colnames(rsid.dt) <- c("SNP1.rsid", "SNP2.rsid")
            chrom.pairs <- cbind(chrom.pairs, rsid.dt)
            chrom.pairs[ , `:=`(raw.score = score)]

        }))

        return(chrom.size.pairs)

    }, results.list = results.list, BPPARAM = bp.param))

    # node weights
    node.weights <- rbindlist(bplapply(seq_along(results.list),
                                       function(d, results.list){

        # chromosome specific results
        chrom.size.res <- results.list[[d]]
        n.top.chroms <- nrow(chrom.size.res)
        chrom.size <- which(colnames(chrom.size.res) == "snp1.rsid") - 1
        these.cols <- seq_len(chrom.size)
        scores <- chrom.size.res$graphical.score

        # snp number data.table
        snp.dt <- chrom.size.res[ , ..these.cols]
        snp.dt$raw.score <- scores

            # snp number data.table
            snp.dt <- chrom.size.res[, ..these.cols]
            snp.dt$raw.score <- scores

            # rsid data.table
            rsid.cols <- seq(chrom.size + 1, 2 * chrom.size)
            rsid.dt <- chrom.size.res[, ..rsid.cols]
            rsid.dt$raw.score <- scores

            # melt
            snp.cols <- seq_len(chrom.size)
            snp.dt.long <- melt(snp.dt, "raw.score", snp.cols)
            snp.dt.long$variable <- NULL
            colnames(snp.dt.long) <- c("raw.score", "SNP")
            rsid.dt.long <- melt(rsid.dt, "raw.score", snp.cols)

            # combine into data.table of SNPs and graphical score
            snp.dt.long$rsid <- rsid.dt.long$value

            # return result
            return(snp.dt.long)
        },
        results.list = results.list, BPPARAM = bp.param
    ))

    # remove the NA's
    edge.weights <- edge.weights[!is.na(edge.weights$raw.score), ]
    node.weights <- node.weights[!is.na(node.weights$raw.score), ]

    # compute score based on score.type
    if (score.type == "max"){

        snp.pair.scores <- edge.weights[ , list(pair.score =
                                                    max(raw.score)),
                                         list(SNP1, SNP2, SNP1.rsid, SNP2.rsid)]
        snp.scores <- node.weights[ , list(snp.score =
                                               max(raw.score)), list(SNP, rsid)]
        setorder(snp.pair.scores, -pair.score)
        setorder(snp.scores, -snp.score)

    } else if (score.type == "sum"){

        snp.pair.scores <- edge.weights[ , list(pair.score =
                                                    sum(raw.score)),
                                         list(SNP1, SNP2, SNP1.rsid, SNP2.rsid)]
        snp.scores <- node.weights[ , list(snp.score = sum(raw.score)),
                                    list(SNP, rsid)]
        setorder(snp.pair.scores, -pair.score)
        setorder(snp.scores, -snp.score)

    } else if (score.type == "logsum"){

        snp.pair.scores <- edge.weights[ , list(pair.score =
                                                    log(sum(raw.score))),
                                         list(SNP1, SNP2, SNP1.rsid, SNP2.rsid)]
        snp.scores <- node.weights[ , list(snp.score = log(sum(raw.score))),
                                    list(SNP, rsid)]
        setorder(snp.pair.scores, -pair.score)
        setorder(snp.scores, -snp.score)
    }

    return(list(pair.scores = snp.pair.scores, snp.scores = snp.scores))
}
