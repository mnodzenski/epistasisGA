#' A function to run a test of the null hypothesis that a collection of SNPs do not exhibit epistasis.
#'
#' This function runs a test of the null hypothesis that a collection of SNPs do not exhibit epistasis.
#'
#' @param snps A character vector of the SNP rsids to be examined (should match up with the rsids in \code{annotation.data}).
#' @param annotation.data A data frame containing columns 'RSID', 'REF' and 'ALT'. Column 'RSID' gives the
#' RSIDs for the input SNPs, with the rows ordered such that the first RSID entry corresponds to the first SNP
#' column in the data passed to function \code{preprocess.genetic.data}, the second RSID corresponds to the second SNP column, etc.
#' @param preprocessed.list The initial list produced by function \code{preprocess.genetic.data}.
#' @param ld.vec A logical matrix indicating whether the the SNPs in \code{snps} should be considered to
#' be in linkage disequilibrium. If not specified, this defaults to the assumption that SNPs located on the
#' same biological chromosome are in LD and those on different chromosomes are not.
#' @param n.different.snps.weight The number by which the number of different SNPs between a case and complement is multiplied in computing the family weights. Defaults to 2.
#' @param n.both.one.weight The number by which the number of SNPs equal to 1 in both the case and complement is multiplied in computing the family weights. Defaults to 1.
#' @param weight.function.int An integer used to assign family weights. Specifically, we use \code{weight.function.int} in a  function that takes the weighted sum
#' of the number of different SNPs and SNPs both equal to one as an argument, denoted as x, and returns a family weight equal to \code{weight.function.int}^x. Defaults to 2.
#' @param n.case.high.risk.thresh The number of cases with the provisional high risk set required to check for recessive patterns of allele inheritance.
#' @param outlier.sd The number of standard deviations from the mean allele count used to determine whether recessive allele coding is appropriate
#' for a given SNP. See the GADGET paper for specific details on the implementation of this argument.
#' @return A list of two elements. Note these two objects will also be written to \code{results.dir}
#' as 'combined.island.results.rds' and 'combined.island.unique.chromosome.results.rds'. Furthermore,
#' this function should not be called twice on the same directory (i.e., only combine the islands one time).
#' \describe{
#'  \item{all.results}{A dataset containing chromosome results across all islands,
#'  where top chromosomes that evolved on multiple distinct islands appear in multiple rows.}
#'  \item{unique.results}{A condensed version of \code{all.results} with one row per distinct chromosome
#'  and an additional variable indicating the number of islands on which that chromosome evolved.}
#' }
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' data(snp.annotations)
#' library(Matrix)
#' chrom.mat <- as.matrix(bdiag(list(matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25))))
#'
#' pp.list <- preprocess.genetic.data(case[, 1:10], father.genetic.data = dad[ , 1:10],
#'                                mother.genetic.data = mom[ , 1:10],
#'                                chrom.mat = chrom.mat[ , 1:10])
#'
#' run.ga(pp.list, n.chromosomes = 4, chromosome.size = 3, results.dir = 'tmp',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'
#' combined.res <- combine.islands('tmp', snp.annotations[ 1:10, ], pp.list)
#'
#' unlink('tmp', recursive = TRUE)
#' unlink('tmp_reg', recursive = TRUE)
#'
#' @importFrom data.table rbindlist setkey setorder `:=`
#' @export

epi.test <- function(snps, annotation.data, preprocessed.list, ld.mat = NULL,
                     n.different.snps.weight = 2, n.both.one.weight = 1,
                     weight.function.int = 2, n.case.high.risk.thresh = 20,
                     outlier.sd = 2.5) {

    ### get the column numbers corresponding to the snps of interest ###
    target <- snps %in% annotation.data$RSID
    if (sum(target) != length(snps)){

        stop("Not all snps are not listed in annotation.data")
    }
    snp.cols <- which(target)
    original.col.numbers <- preprocessed.list$original.col.numbers
    target.snps <- which(original.col.numbers %in% snp.cols)

    ### check for LD, can't do test if all SNPs in ld ###
    if (is.null(ld.mat)){

        ld.mat <- preprocessed.list$chrom.mat[target.snps, target.snps]

    }
    if (all(ld.mat)){

        stop("cannot run test, all SNPs are in LD")

    }

    ### compute the weight lookup table ###
    max.sum <- max(n.different.snps.weight, n.both.one.weight)*length(snps)
    weight.lookup <- vapply(seq_len(max.sum), function(x) weight.function.int^x, 1)

    ### grab pre-processed data ###
    case.genetic.data <- preprocessed.list$case.genetic.data
    complement.genetic.data <- preprocessed.list$complement.genetic.data
    chrom.mat <- preprocessed.list$chrom.mat

    ### Compute matrices of differences between cases and complements ###
    case.minus.comp <- sign(as.matrix(case.genetic.data - complement.genetic.data))
    case.comp.different <- case.minus.comp != 0

    ### Compute matrix indicating whether both the case and control have the same ###
    ### number of copies of the minor allele ###
    both.zero.mat <- complement.genetic.data == 0 & case.genetic.data == 0
    both.one.mat <- complement.genetic.data == 1 & case.genetic.data == 1
    both.two.mat <- complement.genetic.data == 2 & case.genetic.data == 2

    ### compute fitness score ###
    fitness.score <- chrom.fitness.score(case.genetic.data, complement.genetic.data, case.comp.different,
                                    target.snps, case.minus.comp, both.zero.mat, both.one.mat,
                                    both.two.mat, chrom.mat, weight.lookup,
                                    n.different.snps.weight, n.both.one.weight,
                                    n.case.high.risk.thresh, outlier.sd, epi.test = TRUE)

    ### grab the cases and complements with the risk genotype ###
    case.risk <- fitness.score$high.risk.families$cases
    comp.risk <- fitness.score$high.risk.families$complements
    n.case.risk <- length(case.risk)

    ### grab the case data for the risk families ###
    case.risk.data <- case.genetic.data[c(case.risk, comp.risk), target.snps]

    ###  ###


}
