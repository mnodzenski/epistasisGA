#' A function to assign an E-GADGETS (GxGxE) fitness score to a chromosome
#'
#' This function assigns the (currently experimental) E-GADGETS fitness score to a chromosome.
#'
#' @param case.genetic.data The genetic data of the disease affected children
#' from case-parent trios or disease-discordant sibling pairs. If searching for
#' maternal SNPs that are related to risk of disease in the child, some of the
#' columns in \code{case.genetic.data} may contain maternal SNP genotypes
#' (See argument \code{mother.snps} for how to indicate which SNPs columns
#' correspond to maternal genotypes). Columns are SNP allele counts, and rows
#' are individuals. This object may either be of class matrix' OR of class
#' 'big.matrix'. If of class 'big.matrix' it must be file backed as type
#' 'integer' (see the \code{bigmemory} package for more information). The
#' ordering of the columns must be consistent with the LD structure specified
#' in \code{ld.block.vec}. The genotypes cannot be  dosages imputed with
#' uncertainty.
#' @param father.genetic.data The genetic data for the fathers of the cases in
#' \code{case.genetic.data}. This should only be specified when searching for
#' epistasis or GxGxE effects based only on case-parent triads, and not when
#' searching for maternal SNPs that are related to the child's risk of disease.
#' Columns are SNP allele counts, rows are individuals. This object may either
#' be of class 'matrix' OR of class 'big.matrix'. If of class big.matrix' it
#' must be file backed as type 'integer' (see the bigmemory package for more
#' information). The genotypes cannot be dosages imputed with uncertainty.
#' @param mother.genetic.data The genetic data for the mothers of the cases in
#' \code{case.genetic.data}. This should only be specified when searching for
#' epistasis or GxGxE effects based only on case-parent triads, and not when
#' searching for maternal SNPs that are related to the child's risk of disease.
#' Columns are SNP allele counts, rows are individuals. This object may either
#' be of class 'matrix' OR of class 'big.matrix'. If of class big.matrix' it
#' must be file backed as type 'integer' (see the bigmemory package for more
#' information). The genotypes cannot be dosages imputed with uncertainty.
#' @param exposure.mat A matrix of the input categorical
#' and continuous exposures to be used in the experimental
#' E-GADGETS fitness score. If there are categorical exposure variables with
#' more than 2 levels, those should be dummy coded.
#' @param target.snps An integer vector of the columns corresponding to the
#' collection of SNPs, or chromosome, for which the fitness score will be
#' computed.
#' @param weight.lookup A vector that maps a family weight to the weighted sum
#' of the number of different SNPs and SNPs both equal to one.
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
#' @param n.different.snps.weight The number by which the number of different
#' SNPs between a case and complement/unaffected sibling
#' is multiplied in computing the family weights. Defaults to 2.
#' @param n.both.one.weight The number by which the number of SNPs equal to 1 in
#' both the case and complement/unaffected sibling is multiplied in computing
#' the family weights. Defaults to 1.
#' @return A list:
#' \describe{
#'  \item{fitness.score}{The chromosome fitness score.}
#'  \item{sum.dif.vecs}{The element of the Hotelling-Lawley
#'  trace matrix corresponding to each SNP. Larger magnitudes
#'  indicate larger contributions to the score, but are otherwise
#'  difficult to interpret.}
#'  \item{ht_trace}{The Hotelling-Lawley trace statistic from the
#'  transmission-based fitness score component.}
#'  \item{wald_stat}{The Wald statistic from the family-based component
#'  of the fitness score.}
#' }
#'
#' @examples
#'
#' data(case.gxe)
#' data(dad.gxe)
#' data(mom.gxe)
#' data(exposure)
#' case.gxe <- case.gxe + 0.0
#' mom.gxe <- mom.gxe + 0.0
#' dad.gxe <- dad.gxe + 0.0
#' exposure <- as.matrix(exposure + 0.0)
#' weight.lookup <- vapply(seq_len(6), function(x) 2^x, 1.0)
#' res <- GxE.fitness.score(case.gxe, mom.gxe, dad.gxe, exposure, c(1, 4, 7),
#'                           weight.lookup)
#'
#' @export

GxE.fitness.score <- function(case.genetic.data, 
                              mother.genetic.data, 
                              father.genetic.data, 
                              exposure.mat, target.snps, weight.lookup, 
                              null.mean.vec = c(0, 0), null.sd.vec = c(1, 1),
                              n.different.snps.weight = 2, 
                              n.both.one.weight = 1) {

    GxE_fitness_score_mvlm(case.genetic.data, mother.genetic.data, 
                           father.genetic.data, exposure.mat, target.snps, 
                           weight.lookup, null.mean.vec, null.sd.vec, 
                           n.different.snps.weight,
                           n.both.one.weight)

}

