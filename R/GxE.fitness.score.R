#' A function to assign a fitness score to a chromosome
#'
#' This function assigns a fitness score to a chromosome.
#'
#' @param case.genetic.data The genetic data of the disease affected children from case-parent trios or affected/unaffected sibling pairs.
#'  Columns are SNP allele counts, and rows are individuals.
#' The ordering of the columns must be consistent with the LD structure specified in \code{ld.block.vec}.
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
#' @param ld.block.vec An integer vector specifying the linkage blocks of the input SNPs. As an example, for 100 candidate SNPs, suppose
#' we specify \code{ld.block.vec <- c(25, 75, 100)}. This vector indicates that the input genetic data has 3 distinct linkage blocks, with
#' SNPs 1-25 in the first linkage block, 26-75 in the second block, and 76-100 in the third block. Note that this means the ordering of the columns (SNPs)
#' in \code{case.genetic.data} must be consistent with the LD blocks specified in \code{ld.block.vec}. In the absence of outside information,
#' a reasonable default is to consider SNPs to be in LD if they are located on the same biological chromosome. If not specified, this defaults
#' to assuming all input SNPs are in linkage, which may be overly conservative and could adversely affect performance.
#' @param weight.lookup A vector that maps a family weight to the weighted sum of the number of different SNPs and SNPs both equal to one.
#' @param case2.mat A logical matrix indicating whether, for each SNP, the case carries 2 copies of the minor allele.
#' @param case0.mat A logical matrix indicating whether, for each SNP, the case carries 0 copies of the minor allele.
#' @param n.different.snps.weight The number by which the number of different SNPs between a case and complement/unaffected sibling
#'  is multiplied in computing the family weights. Defaults to 2.
#' @param n.both.one.weight The number by which the number of SNPs equal to 1 in both the case and complement/unaffected sibling
#'  is multiplied in computing the family weights. Defaults to 1.
#' @param recode.threshold For a given SNP, the minimum test statistic required to recode and recompute the fitness score using recessive coding. Defaults to 3.
#' See the GADGETS paper for specific details.
#' @param exposure A categorical factor vector corresponding to environmental exposures of the cases.
#' @return A list:
#' \describe{
#'  \item{fitness.score}{The chromosome fitness score.}
#'  \item{sum.dif.vecs}{The weighted mean difference vector corresponding to the chromosome,
#'  with each element divided by it's pseudo-standard error.
#'  The magnitudes of these values are not particularly important, but the sign is useful.
#'  A positive value for a given SNP indicates the minor allele is positively associated with
#'  disease status, while a negative value implies the reference (‘wild type’) allele is
#'  positively associated with the disease.}
#'  \item{q}{The fraction of provisional risk alleles carried by cases with the full risk set
#'  over the total number of risk alleles carried by either a case or complement with the full risk set.}
#'  \item{pseudo.t2}{The pseudo T^2 value for the chromosome.}
#'  \item{risk.set.alleles}{A vector indicating the number risk alleles a case or complement must have
#'   for each SNP in \code{target.snps} for the case or complement to be classified as having the
#'   proposed risk set. '1+' indicates at least one copy of the risk allele is required, while '2'
#'   indicates 2 copies are needed. The risk allele can be determined based on the signs of the elements
#'   of \code{sum.dif.vecs}, where a negative value indicates the major allele for a given SNP is
#'   the risk allele, while a positive value implicates the minor allele.}
#'   \item{inf.families}{An integer vector of the informative family rows. Only returned if \code{epi.test} = TRUE.}
#'   \item{xbar}{The weighted mean difference sign vector. Only returned if \code{GxE} is TRUE.}
#'   \item{sigma.hat}{The pseudo covariance matrix of the weighted mean difference sign vector. Only returned if \code{GxE} is TRUE.}
#'   \item{w}{The sum of the family weights. Only returned if \code{GxE} is TRUE.}
#' }
#'
#' @examples
#'
#' data(case.gxe)
#' data(dad.gxe)
#' data(mom.gxe)
#' data(exposure)
#' exposure <- as.integer(exposure)
#' comp.gxe <- mom.gxe + dad.gxe - case.gxe
#' case.comp.diff <- case.gxe != comp.gxe
#' case.minus.comp <- case.gxe - comp.gxe
#' both.one.mat <- case.gxe == 1 & comp.gxe == 1
#' case2.mat <- case.gxe == 2
#' case0.mat <- case.gxe == 0
#' comp2.mat <- comp.gxe == 2
#' comp0.mat <- comp.gxe == 0
#' ld.block.vec <- cumsum(rep(25, 4))
#' weight.lookup <- vapply(seq_len(6), function(x) 2^x, 1)
#' case.list <- split(case.gxe, exposure)
#' exposure.levels <- as.integer(names(case.list))
#' case.list <- lapply(case.list, as.matrix)
#' comp.list <- split(comp.gxe, exposure)
#' comp.list <- lapply(comp.list, as.matrix)
#' case.comp.diff.list <- split(case.comp.diff, exposure)
#' case.comp.diff.list <- lapply(case.comp.diff.list, as.matrix)
#' case.minus.comp.list <- split(case.minus.comp, exposure)
#' case.minus.comp.list <- lapply(case.minus.comp.list, as.matrix)
#' both.one.mat.list <- split(both.one.mat, exposure)
#' both.one.mat.list <- lapply(both.one.mat.list, as.matrix)
#' case2.mat.list <- split(case2.mat, exposure)
#' case2.mat.list <- lapply(case2.mat.list, as.matrix)
#' case0.mat.list <- split(case0.mat, exposure)
#' case0.mat.list <- lapply(case0.mat.list, as.matrix)
#' comp2.mat.list <- split(comp2.mat, exposure)
#' comp2.mat.list <- lapply(comp2.mat.list, as.matrix)
#' comp0.mat.list <- split(comp0.mat, exposure)
#' comp0.mat.list <- lapply(comp0.mat.list, as.matrix)
#'
#' GxE.fitness.score(case, comp, case.comp.diff, c(1, 4, 7),
#'                     case.minus.comp, both.one.mat,
#'                     ld.block.vec, weight.lookup,
#'                     case2.mat, case0.mat, exposure)
#'
#' @export

GxE.fitness.score <- function(case.genetic.data, complement.genetic.data, case.comp.differences,
                                target.snps, cases.minus.complements, both.one.mat,
                                ld.block.vec, weight.lookup, case2.mat, case0.mat, exposure,
                                n.different.snps.weight = 2, n.both.one.weight = 1,
                                recode.threshold = 3) {



}

