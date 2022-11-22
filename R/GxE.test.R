#' A function to run a test of interaction between a collection of SNPs and exposures (experimental).
#'
#' This function runs a permutation based test run a test of
#' interaction between a collection of SNPs and exposure variables.
#'
#' @param snp.cols An integer vector specifying the columns in the input data
#' containing the SNPs to be tested.
#' @param preprocessed.list The initial list produced by function
#' \code{preprocess.genetic.data}.
#' @param null.mean.vec A vector of estimated null means for each
#' of the components of the E-GADGETS fitness score.
#' It should be set to the values of the "null.mean" element of the file
#' "null.mean.sd.info.rds" for the observed data, that is saved by the
#' \code{run.gadgets} function.
#' @param null.sd.vec A vector of estimated null means for each
#' of the components of the E-GADGETS fitness score.
#' It should be set to the values of the "null.se" element of the file
#' "null.mean.sd.info.rds" for the observed data, that is saved by the
#' \code{run.gadgets} function.
#' @param n.permutes The number of permutations on which to base the test.
#' Defaults to 10000.
#' @param n.different.snps.weight The number by which the number of different
#' SNPs between a case and complement/unaffected sibling is multiplied in
#' computing the family weights. Defaults to 2.
#' @param n.both.one.weight The number by which the number of SNPs equal to 1
#' in both the case and complement/unaffected sibling is multiplied in computing
#' the family weights. Defaults to 1.
#' @param weight.function.int An integer used to assign family weights.
#' Specifically, we use \code{weight.function.int} in a function that takes the
#' weighted sum of the number of different SNPs and SNPs both equal to one as an
#' argument, denoted as x, and returns a family weight equal to
#' \code{weight.function.int}^x. Defaults to 2.
#' @return A list of three elements:
#' \describe{
#'  \item{pval}{The p-value of the test.}
#'  \item{obs.fitness.score}{The fitness score from the observed data}
#'  \item{perm.fitness.scores}{A vector of fitness scores for the
#'   permuted datasets.}
#' }
#' @examples
#'
#' data(case.gxe)
#' data(dad.gxe)
#' data(mom.gxe)
#' data(exposure)
#' data(snp.annotations.mci)
#' pp.list <- preprocess.genetic.data(case.gxe, father.genetic.data = dad.gxe,
#'                                mother.genetic.data = mom.gxe,
#'                                ld.block.vec = rep(6, 4),
#'                                categorical.exposures = exposure)
#'  run.gadgets(pp.list, n.chromosomes = 5, chromosome.size = 3,
#'              results.dir = "tmp_gxe", cluster.type = "interactive",
#'              registryargs = list(file.dir = "tmp_reg_gxe", seed = 1300),
#'              n.islands = 8, island.cluster.size = 4,
#'              n.migrations = 1)
#'
#'
#' combined.res <- combine.islands('tmp_gxe', snp.annotations.mci, pp.list, 1)
#' top.snps <- as.vector(t(combined.res[1, 1:3]))
#' set.seed(10)
#' GxE.test.res <- GxE.test(top.snps, pp.list)
#'
#' unlink('tmp_gxe', recursive = TRUE)
#' unlink('tmp_reg_gxe', recursive = TRUE)
#'
#' @export

GxE.test <- function(snp.cols, preprocessed.list, null.mean.vec = c(0, 0),
                     null.sd.vec =c(1,1), n.permutes = 10000,
                     n.different.snps.weight = 2, n.both.one.weight = 1,
                     weight.function.int = 2) {

    # run the test via cpp
    GxE_test(snp.cols, preprocessed.list, null.mean.vec, null.sd.vec, n.permutes,
             n.different.snps.weight, n.both.one.weight, weight.function.int)

}
