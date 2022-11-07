#' A function to run a test of the null hypothesis that a collection of SNPs do not exhibit epistasis, conditional
#' upon observed marginal SNP-disease associations.
#'
#' This function runs a permutation based test of the null hypothesis that a collection of SNPs do not exhibit epistasis,
#' conditional upon observed marginal SNP-disease associations.
#'
#' @param snp.cols An integer vector specifying the columns in the input data containing the SNPs to be tested.
#' @param preprocessed.list The initial list produced by function \code{preprocess.genetic.data}.
#' @param n.permutes The number of permutations on which to base the test. Defaults to 1000.
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
#' @param maternal.fetal.test A boolean indicating whether the test specifically for a maternal-fetal interaction should be run. Defaults to FALSE.
#' @return A list of thee elements:
#' \describe{
#'  \item{pval}{The p-value of the test.}
#'  \item{obs.fitness.score}{The fitness score from the observed data}
#'  \item{perm.fitness.scores}{A vector of fitness scores for the permuted datasets.}
#' }
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' data(snp.annotations)
#'
#' pp.list <- preprocess.genetic.data(as.matrix(case),
#'                                father.genetic.data = as.matrix(dad),
#'                                mother.genetic.data = as.matrix(mom),
#'                                ld.block.vec = rep(25, 4))
#'
#' run.gadgets(pp.list, n.chromosomes = 5, chromosome.size = 3,
#'        results.dir = "tmp", cluster.type = "interactive",
#'        registryargs = list(file.dir = "tmp_reg", seed = 1300),
#'        n.islands = 8, island.cluster.size = 4,
#'        n.migrations = 2)
#'
#' combined.res <- combine.islands('tmp', snp.annotations, pp.list, 2)
#'
#' top.snps <- as.vector(t(combined.res[1, 1:3]))
#' set.seed(10)
#' epi.test.res <- epistasis.test(top.snps, pp.list)
#'
#' unlink('tmp', recursive = TRUE)
#' unlink('tmp_reg', recursive = TRUE)
#' @export

epistasis.test <- function(snp.cols, preprocessed.list, n.permutes = 10000,
                     n.different.snps.weight = 2, n.both.one.weight = 1,
                     weight.function.int = 2, recessive.ref.prop = 0.75,
                     recode.test.stat = 1.64, maternal.fetal.test = FALSE) {

    # run the epistasis test via cpp
    epistasis_test(snp.cols, preprocessed.list, n.permutes,
                   n.different.snps.weight, n.both.one.weight, weight.function.int,
                   recessive.ref.prop, recode.test.stat, TRUE, maternal.fetal.test)

}
