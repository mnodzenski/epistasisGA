#' A function to assign a fitness score to a chromosome
#'
#' This function assigns a fitness score to a chromosome. It is a wrapper for the Rcpp function chrom_fitness_score.
#'
#' @param case.genetic.data The genetic data of the disease affected children from case-parent trios or affected/unaffected sibling pairs.
#'  Columns are SNP allele counts, and rows are individuals.
#' The ordering of the columns must be consistent with the LD structure specified in \code{block.ld.mat}.
#' @param complement.genetic.data A genetic dataset from the complements of the cases, where
#' \code{complement.genetic.data} = mother SNP counts + father SNP counts - case SNP counts.
#' Columns are SNP allele counts, rows are families. If using affected/unaffected sibling pairs, this should contain
#' the unaffected sibling genotypes.
#' @param case.comp.differences A data frame or matrix indicating \code{case.genetic.data} != \code{complement.genetic.data},
#' where rows correspond to individuals and columns correspond to snps.
#' @param target.snps A numeric vector of the columns corresponding to the collection of SNPs, or chromosome, for which the fitness score will be computed.
#' @param cases.minus.complements A matrix equal to \code{case.genetic.data} - \code{complement genetic data}.
#' @param both.one.mat A matrix whose elements indicate whether both the case and complement have one copy of the minor allele,
#' equal to \code{case.genetic.data == 1 & complement.genetic.data == 1}.
#' @param block.ld.mat A logical, block diagonal matrix indicating whether the SNPs in \code{case.genetic.data} should be considered
#'  to be in linkage disequilibrium. Note that this means the ordering of the columns (SNPs) in \code{case.genetic.data} must be consistent
#'  with the LD blocks specified in \code{ld.block.mat}. In the absence of outside information, a reasonable default is to consider SNPs
#'  to be in LD if they are located on the same biological chromosome.
#' @param weight.lookup A vector that maps a family weight to the weighted sum of the number of different SNPs and SNPs both equal to one.
#' @param case2.mat A logical matrix indicating whether, for each SNP, the case carries 2 copies of the minor allele.
#' @param case0.mat A logical matrix indicating whether, for each SNP, the case carries 0 copies of the minor allele.
#' @param comp2.mat A logical matrix indicating whether, for each SNP, the complement/unaffected sibling carries 2 copies of the minor allele.
#' @param comp0.mat A logical matrix indicating whether, for each SNP, the complement/unaffected sibling carries 0 copies of the minor allele.
#' @param n.different.snps.weight The number by which the number of different SNPs between a case and complement/unaffected sibling
#'  is multiplied in computing the family weights. Defaults to 2.
#' @param n.both.one.weight The number by which the number of SNPs equal to 1 in both the case and complement/unaffected sibling
#'  is multiplied in computing the family weights. Defaults to 1.
#' @param recessive.ref.prop The proportion to which the observed proportion of informative cases with the provisional risk genotype(s) will be compared
#' to determine whether to recode the SNP as recessive. Defaults to 0.75.
#' @param recode.test.stat For a given SNP, the minimum test statistic required to recode and recompute the fitness score using recessive coding. Defaults to 1.64.
#' See the GADGETS paper for specific details.
#' @param epi.test A logical indicating whether the function should return the information required to run function \code{epistasis.test}.
#' for a given SNP. See the GADGETS paper for specific details on the implementation of this argument.
#' @param dif.coding A logical indicating whether, for a given SNP, the case - complement genotype difference should
#' be coded as the sign of the difference (defaulting to false) or the raw difference.
#' @param GxE a logical indicating whether the function is being used to compute a GxE fitness score, defaulting to false.
#' @return A list:
#' \describe{
#'  \item{fitness_score}{The chromosome fitness score.}
#'  \item{sum_dif_vecs}{The weighted mean difference vector corresponding to the chromosome,
#'  with each element divided by it's pseudo-standard error.
#'  The magnitudes of these values are not particularly important, but the sign is useful.
#'  A positive value for a given SNP indicates the minor allele is positively associated with
#'  disease status, while a negative value implies the reference (‘wild type’) allele is
#'  positively associated with the disease.}
#'  \item{q}{The fraction of provisional risk alleles carried by cases with the full risk set
#'  over the total number of risk alleles carried by either a case or complement with the full risk set.}
#'  \item{risk_set_alleles}{A vector indicating the number risk alleles a case or complement must have
#'   for each SNP in \code{target.snps} for the case or complement to be classified as having the
#'   proposed risk set. '1+' indicates at least one copy of the risk allele is required, while '2'
#'   indicates 2 copies are needed. The risk allele can be determined based on the signs of the elements
#'   of \code{sum_dif_vecs}, where a negative value indicates the major allele for a given SNP is
#'   the risk allele, while a positive value implicates the minor allele.}
#'  \item{inf_families}{An integer vector of the informative family rows. Only returned if \code{epi.test} = TRUE.}
#' }
#'
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' case <- as.matrix(case)
#' dad <- as.matrix(dad)
#' mom <- as.matrix(mom)
#' comp <- mom + dad - case
#' case.comp.diff <- case != comp
#' case.minus.comp <- sign(case - comp)
#' storage.mode(case.minus.comp) <- "integer"
#' both.one.mat <- case == 1 & comp == 1
#' case2.mat <- case == 2
#' case0.mat <- case == 0
#' comp2.mat <- comp == 2
#' comp0.mat <- comp == 0
#' library(Matrix)
#' block.ld.mat <- as.matrix(bdiag(list(matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25))))
#' weight.lookup <- vapply(seq_len(6), function(x) 2^x, 1)
#' storage.mode(weight.lookup) <- "integer"
#' chrom.fitness.score(case, comp, case.comp.diff, c(1, 4, 7),
#'                     case.minus.comp, both.one.mat,
#'                     block.ld.mat, weight.lookup,
#'                     case2.mat, case0.mat, comp2.mat,
#'                     comp0.mat)
#'
#' @export

chrom.fitness.score <- function(case.genetic.data, complement.genetic.data, case.comp.differences,
                                target.snps, cases.minus.complements, both.one.mat,
                                block.ld.mat, weight.lookup, case2.mat, case0.mat,
                                comp2.mat, comp0.mat, n.different.snps.weight = 2,
                                n.both.one.weight = 1, recessive.ref.prop = 0.75,
                                recode.test.stat = 1.64, epi.test = FALSE, dif.coding = FALSE,
                                GxE = FALSE) {

  chrom_fitness_score(case.genetic.data, complement.genetic.data, case.comp.differences, target.snps,
                      cases.minus.complements, both.one.mat, block.ld.mat, weight.lookup, case2.mat,
                      case0.mat, comp2.mat, comp0.mat, n.different.snps.weight, n.both.one.weight,
                      recessive.ref.prop, recode.test.stat, epi.test, dif.coding, GxE)

}

