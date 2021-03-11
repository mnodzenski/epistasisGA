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
#' @param pval.thresh A numeric value between 0 and 1 specifying the epistasis test p-value threshold for a
#' chromosome to contribute to the network. Any chromosomes with epistasis p-value greater than \code{pval.thresh}
#' will not contribute to network plots. The argument defaults to 0.05.
#' @param n.permutes The number of permutations on which to base the epistasis tests. Defaults to 10000.
#' @param n.different.snps.weight The number by which the number of different SNPs between a case and complement/unaffected sibling
#'  is multiplied in computing the family weights. Defaults to 2.
#' @param n.both.one.weight The number by which the number of SNPs equal to 1 in both the case and complement/unaffected sibling
#' is multiplied in computing the family weights. Defaults to 1.
#' @param weight.function.int An integer used to assign family weights. Specifically, we use \code{weight.function.int} in a function that takes the weighted sum
#' of the number of different SNPs and SNPs both equal to one as an argument, denoted as x, and
#' returns a family weight equal to \code{weight.function.int}^x. Defaults to 2.
#' @param recessive.ref.prop The proportion to which the observed proportion of informative cases with the provisional risk genotype(s) will be compared
#' to determine whether to recode the SNP as recessive. Defaults to 0.75.
#' @param recode.test.stat For a given SNP, the minimum test statistic required to recode and recompute the fitness score using recessive coding. Defaults to 1.64.
#' See the GADGETS paper for specific details.
#' @param dif.coding A logical indicating whether, for a given SNP, the case - complement genotype difference should
#' be coded as the sign of the difference (defaulting to false) or the raw difference.
#' @param bp.param The BPPARAM argument to be passed to bplapply. See \code{BiocParallel::bplapply} for more details.
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
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  combined.res2 <- combine.islands('tmp_2', snp.annotations[ target.snps, ], pp.list, 2)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#'  #observed data chromosome size 3
#'  run.gadgets(pp.list, n.chromosomes = 5, chromosome.size = 3, results.dir = 'tmp_3',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'  combined.res3 <- combine.islands('tmp_3', snp.annotations[ target.snps, ], pp.list, 2)
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
#' @importFrom BiocParallel bplapply bpparam
#' @export

compute.pair.scores <- function(results.list, pp.list, n.top.chroms = 50, score.type = "logsum",
                                pval.thresh = 0.05, n.permutes = 10000, n.different.snps.weight = 2,
                                n.both.one.weight = 1, weight.function.int = 2, recessive.ref.prop = 0.75,
                                recode.test.stat = 1.64, dif.coding = FALSE, bp.param = bpparam()) {

    ## make sure we have the correct number of chromosomes in each element of the results list
    ## and then return just the chromosomes of interest
    chrom.list <- lapply(results.list, function(chrom.size.data){


        n.obs.chroms <- sum(!is.na(chrom.size.data$fitness.score))
        if (n.obs.chroms != n.top.chroms){

            stop("Each element of results.list must contain n.top.scores fitness scores.")

        }
        chrom.size <- sum(grepl("snp", colnames(chrom.size.data)))/5
        these.cols <- seq_len(chrom.size)
        chrom.mat <- as.matrix(chrom.size.data[ , ..these.cols])
        chrom.list <- split(chrom.mat, seq_len(nrow(chrom.mat)))
        return(chrom.list)

    })

    ## compute graphical scores based on epistasis test
    n2log.epi.pvals <- bplapply(chrom.list, function(chrom.size.list, pp.list, n.permutes,
                                                       n.different.snps.weight, n.both.one.weight,
                                                       weight.function.int, recessive.ref.prop,
                                                       recode.test.stat, dif.coding){

        n2log_epistasis_pvals(chrom.size.list, pp.list, n.permutes,
                              n.different.snps.weight, n.both.one.weight, weight.function.int,
                              recessive.ref.prop, recode.test.stat, dif.coding)

    }, pp.list = pp.list, n.permutes = n.permutes, n.different.snps.weight = n.different.snps.weight,
    n.both.one.weight = n.both.one.weight, weight.function.int = weight.function.int,
    recessive.ref.prop = recessive.ref.prop, recode.test.stat = recode.test.stat,
    dif.coding = dif.coding, BPPARAM = bp.param)

   ## add those scores to the obs data
   for (i in seq_along(n2log.epi.pvals)){

       results.list[[i]]$graphical.score <- n2log.epi.pvals[[i]]
       results.list[[i]] <- results.list[[i]][results.list[[i]]$graphical.score >= -2*log(pval.thresh), ]

   }

    #make sure we have some edges
    n.edges <- vapply(results.list, nrow, 1)
    if (sum(n.edges) == 0){
        stop("No SNP pairs meet p-value threshold")
    }

    all.edge.weights <- rbindlist(bplapply(seq_along(results.list), function(d, results.list){

        chrom.size.res <- results.list[[d]]
        chrom.size <- sum(grepl("snp", colnames(chrom.size.res)))/5
        rbindlist(lapply(seq_len(n.top.chroms), function(res.row) {

            chrom.res <- chrom.size.res[res.row, ]
            score <- chrom.res$graphical.score
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
            chrom.pairs[ , `:=`(raw.score = score)]

            #take average score for each pair
            out.dt <- chrom.pairs[ , list(graphical.score = mean(score)),
                                   list(SNP1, SNP2, SNP1.rsid, SNP2.rsid)]
            return(out.dt)

        }))

    }, results.list = results.list, BPPARAM = bp.param))

    # remove the NA's
    all.edge.weights <- all.edge.weights[!is.na(all.edge.weights$graphical.score), ]

    # compute edge score based on score.type
    if (score.type == "max"){

        out.dt <- all.edge.weights[ , list(edge.score = max(graphical.score)), list(SNP1, SNP2, SNP1.rsid, SNP2.rsid)]
        setorder(out.dt, -edge.score)

    } else if (score.type == "sum"){

        out.dt <- all.edge.weights[ , list(edge.score = sum(graphical.score)),
                                    list(SNP1, SNP2, SNP1.rsid, SNP2.rsid)]
        setorder(out.dt, -edge.score)

    } else if (score.type == "logsum"){

        out.dt <- all.edge.weights[ , list(edge.score = log(1 + (sum(graphical.score)))),
                                    list(SNP1, SNP2, SNP1.rsid, SNP2.rsid)]
        setorder(out.dt, -edge.score)

    }

    return(out.dt)

}


