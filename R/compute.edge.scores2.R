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
#' @return A data.table where the first four columns represent SNPs and the fifth column (edge.score)
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
#'                                block.ld.mat = block.ld.mat[1:10 , 1:10])
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
#'  ## create list of results
#'  final.results <- list(combined.res2$unique.results[1:3, ], combined.res3$unique.results[1:3, ])
#'
#'  ## compute edge scores
#'  edge.dt <- compute.edge.scores2(final.results, pp.list, 3)
#'
#'  lapply(c('tmp_2', 'tmp_3'), unlink, recursive = TRUE)
#'
#' @importFrom data.table data.table rbindlist setorder
#' @importFrom matrixStats colSds
#' @importFrom utils combn
#' @importFrom BiocParallel bplapply SerialParam
#' @export

compute.edge.scores2 <- function(results.list, pp.list, n.top.chroms = 50, score.type = "logsum",
                                 bp.param = SerialParam(), zscore.thresh = NULL) {

    ## make sure we have the correct number of chromosomes in each element of the results list
    n.chroms.vec <- lapply(results.list, function(chrom.size.list){

        n.obs.chroms <- sum(!is.na(chrom.size.list$fitness.score))
        if (n.obs.chroms != n.top.chroms){

            stop("Each element of results.list must contain n.top.scores fitness scores.")

        }

    })

    # find the correct columns in the pre-processed data
    case <- pp.list$case.genetic.data
    comp <- pp.list$complement.genetic.data

    ## compute re-scaled fitness score based on ranks compared to permutes
    all.edge.weights <- rbindlist(bplapply(results.list, function(obs.res){

        chrom.size <- sum(grepl("snp", colnames(obs.res)))/5

        # pick out chromosome col numbers
        chrom.num.cols <- seq_len(chrom.size)
        chrom.nums <- obs.res[, ..chrom.num.cols]

        # pick out the diff vecs
        dif.vec.cols <- seq(chrom.size*3 + 1, chrom.size*4)
        dif.vecs <- obs.res[, ..dif.vec.cols]

        # pick out rsids
        rsid.cols <- seq(chrom.size + 1, chrom.size*2)
        rsids <- obs.res[, ..rsid.cols]

        # pick out allele copies
        allele.copy.cols <- seq(chrom.size*4 + 1, chrom.size*5)
        allele.copies <- obs.res[, ..allele.copy.cols]

        # loop over chromosomes
        pair.idx.mat <- t(combn(seq_len(chrom.size), 2))
        scored.pairs <- rbindlist(lapply(seq_len(n.top.chroms), function(x){

            target.chrom <- as.vector(t(chrom.nums[x, ]))
            target.cols  <- match(target.chrom, pp.list$original.col.numbers)
            target.chrom.mat <- t(apply(pair.idx.mat, 1, function(y) target.chrom[y]))

            target.dv <- as.vector(t(dif.vecs[x, ]))
            target.dv.mat <- t(apply(pair.idx.mat, 1, function(y) target.dv[y]))

            target.rsid <- as.vector(t(rsids[x, ]))
            target.rsid.mat <- t(apply(pair.idx.mat, 1, function(y) target.rsid[y]))

            target.ac <- as.vector(t(allele.copies[x, ]))
            target.ac.mat <- t(apply(pair.idx.mat, 1, function(y) target.ac[y]))

            #figure out who has all risk genotypes
            risk.genotypes <- ifelse(target.dv <= 0 & target.ac == "1+", 1,
                                     ifelse(target.dv <= 0 & target.ac == "2", 0,
                                            ifelse(target.dv > 0 & target.ac == "1+", 1, 2)))

            risk.dirs <- sign(target.dv)
            pos.risk <- which(risk.dirs > 0)
            pos.risk.cols <- target.cols[pos.risk]
            print(pos.risk.cols)
            neg.risk <- which(risk.dirs <= 0)
            neg.risk.cols <- target.cols[neg.risk]
            n.pos <- length(pos.risk)
            n.neg <- length(neg.risk)

            if (n.neg > 0) {

                neg.risk.genotypes <- risk.genotypes[neg.risk]
                neg.risk.geno.mat <- matrix(rep(neg.risk.genotypes, each = nrow(case)),
                                            ncol = n.neg, byrow = FALSE)
                n1 <- rowSums(case[, neg.risk.cols, drop = FALSE] <= neg.risk.geno.mat)
                n2 <- rowSums(comp[, neg.risk.cols, drop = FALSE] <= neg.risk.geno.mat)
            }
            if (n.pos > 0) {

                pos.risk.genotypes <- risk.genotypes[pos.risk]
                pos.risk.geno.mat <- matrix(rep(pos.risk.genotypes, each = nrow(case)),
                                            ncol = n.pos, byrow = FALSE)
                p1 <- rowSums(case[, pos.risk.cols, drop = FALSE] <= pos.risk.geno.mat)
                p2 <- rowSums(comp[, pos.risk.cols, drop = FALSE] <= pos.risk.geno.mat)

            }
            if (n.pos == 0) {
                case.high.risk <- n1 == chrom.size
                comp.high.risk <- n2 == chrom.size
            } else if (n.neg == 0) {
                case.high.risk <- p1 == chrom.size
                comp.high.risk <- p2 == chrom.size
            } else {
                case.high.risk <- p1 + n1 == chrom.size
                comp.high.risk <- p2 + n2 == chrom.size
            }

            both.high.risk <- case.high.risk & comp.high.risk
            these.families <- (case.high.risk | comp.high.risk) & !both.high.risk
            n.families <- sum(these.families)

            # looping over pairs
            if ( n.families > 0){

                pair.score <- vapply(seq_len(nrow(target.chrom.mat)), function(z){

                    # pick out the chromosome
                    target.pair <- target.chrom.mat[z, ]
                    target.pair.dv <- target.dv.mat[z, ]
                    target.pair.ac <- target.ac.mat[z, ]

                    # determine risk genotypes
                    risk.genotypes <- ifelse(target.pair.dv <= 0 & target.pair.ac == "1+", 1,
                                             ifelse(target.pair.dv <= 0 & target.pair.ac == "2", 0,
                                                    ifelse(target.pair.dv > 0 & target.pair.ac == "1+", 1, 2)))

                    target.pair.cols <- match(target.pair, pp.list$original.col.numbers)

                    # determine cases and comps with risk genotypes
                    risk.dirs <- sign(target.pair.dv)
                    pos.risk <- which(risk.dirs > 0)
                    pos.risk.cols <- target.pair.cols[pos.risk]
                    neg.risk <- which(risk.dirs <= 0)
                    neg.risk.cols <- target.pair.cols[neg.risk]
                    n.pos <- length(pos.risk)
                    n.neg <- length(neg.risk)

                    if (n.neg > 0) {

                        neg.risk.genotypes <- risk.genotypes[neg.risk]
                        neg.risk.geno.mat <- matrix(rep(neg.risk.genotypes, each = n.families),
                                                    ncol = n.neg, byrow = FALSE)
                        n1 <- rowSums(case[these.families, neg.risk.cols, drop = FALSE] <= neg.risk.geno.mat)
                        n2 <- rowSums(comp[these.families, neg.risk.cols, drop = FALSE] <= neg.risk.geno.mat)
                    }
                    if (n.pos > 0) {

                        pos.risk.genotypes <- risk.genotypes[pos.risk]
                        pos.risk.geno.mat <- matrix(rep(pos.risk.genotypes, each = n.families),
                                                    ncol = n.pos, byrow = FALSE)
                        p1 <- rowSums(case[these.families, pos.risk.cols, drop = FALSE] <= pos.risk.geno.mat)
                        p2 <- rowSums(comp[these.families, pos.risk.cols, drop = FALSE] <= pos.risk.geno.mat)

                    }
                    if (n.pos == 0) {
                        case.high.risk <- n1 == 2
                        comp.high.risk <- n2 == 2
                    } else if (n.neg == 0) {
                        case.high.risk <- p1 == 2
                        comp.high.risk <- p2 == 2
                    } else {
                        case.high.risk <- p1 + n1 == 2
                        comp.high.risk <- p2 + n2 == 2
                    }

                    both.high.risk <- case.high.risk & comp.high.risk
                    n.case.high.risk <- sum(case.high.risk[!both.high.risk])
                    n.comp.high.risk <- sum(comp.high.risk[!both.high.risk])
                    n <- n.case.high.risk + n.comp.high.risk
                    phat <- n.case.high.risk/n
                    if (n == 0){
                        phat <- 0.5
                    } else if (phat == 1){
                        phat <- 1-(10^-8)
                    } else if (phat == 0){
                        phat <- 10^-8
                    }
                    zscore <- (phat - 0.5)/sqrt(phat*(1-phat)/n)
                    return(zscore)

                }, 1.0)

            } else {

                pair.score <- 0.5

            }

            target.pair.dt <- data.table(cbind(target.chrom.mat, target.rsid.mat))
            target.pair.dt$score <- pair.score
            if (!is.null(zscore.thresh)){

                target.pair.dt <- target.pair.dt[target.pair.dt$score >= zscore.thresh, ]

            }

            return(target.pair.dt)

        }))
        colnames(scored.pairs) <- c("SNP1", "SNP2", "SNP1.rsid", "SNP2.rsid", "score")

        #average over pairs
        final.scores <- scored.pairs[ , list(pair.score = mean(score)),
                                     list(SNP1, SNP2, SNP1.rsid, SNP2.rsid)]
        return(final.scores)

    }, BPPARAM = bp.param))

    # compute edge score based on score.type
    if (score.type == "max"){

        out.dt <- all.edge.weights[ , list(edge.score = max(pair.score)), list(SNP1, SNP2, SNP1.rsid, SNP2.rsid)]
        setorder(out.dt, -edge.score)

    } else if (score.type == "sum"){

        out.dt <- all.edge.weights[ , list(edge.score = sum(pair.score)),
                                    list(SNP1, SNP2, SNP1.rsid, SNP2.rsid)]
        out.dt <- out.dt[out.dt$edge.score > 0, ]

        setorder(out.dt, -edge.score)

    } else if (score.type == "logsum"){

        sum.dt <- all.edge.weights[ , list(sum.score = sum(pair.score)),
                                    list(SNP1, SNP2, SNP1.rsid, SNP2.rsid)]
        sum.dt <- sum.dt[sum.dt$sum.score > 0, ]
        out.dt <- sum.dt[ , list(edge.score = log(1 + sum.score)),
                          list(SNP1, SNP2, SNP1.rsid, SNP2.rsid)]
        setorder(out.dt, -edge.score)

    }

    return(out.dt)

}


