#' A function to generate a fitness score for a subset of snps from a dataset of case genetic markers
#'
#' This function returns a fitness score for a set of snps, where the fitness score is the squared vector length of the weighted sum of difference vectors between cases and complements.
#'
#' @param case.comp.differences a data frame or matrix indicating case genetic data != complement genetic data, where rows correspond to individuals and columns correspond to snps.
#' @param target.snps A numeric vector of the columns corresponding to the snps for which the fitness score will be computed.
#' @param cases.minus.complements A matrix equal to case genetic data - complement genetic data.
#' @param both.one.mat A matrix whose elements indicate whether both the case and complement have one copy of the alternate allele, equal to (case.genetic.data == 1 & complement.genetic.data == 1).
#' @param n.different.snps.weight The number by which the number different snps between case and control is multiplied in computing the family weights. Defaults to 2.
#' @param n.both.one.weight The number by which the number of snps equal to 1 in both case and control is multiplied in computing the family weights. Defaults to 1.
#' @return A scalar, the fitness score for the given set of snps.
#'
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' comp <- mom + dad - case
#' case.comp.diff <- case != comp
#' case.minus.comp <- case - comp
#' both.one.mat <- case == 1 & comp == 1
#' chrom.fitness.score(case.comp.diff, target.snps = c(1, 4, 7), case.minus.comp, both.one.mat)
#'
#' @export

chrom.fitness.score <- function(case.comp.differences, target.snps, cases.minus.complements, both.one.mat,
                                n.different.snps.weight = 2, n.both.one.weight = 1){

  ### pick out the differences for the target snps ###
  case.comp.diff <- case.comp.differences[ , target.snps]

  ### determine whether families are informative for the set of target.snps ###
  total.different.snps <-rowSums(case.comp.diff)
  informative.families <- total.different.snps != 0
  n.informative.families <- sum(informative.families)

  ### compute weights ###
  both.one <- rowSums(both.one.mat[informative.families , target.snps])
  family.weights <- n.both.one.weight*both.one + n.different.snps.weight*total.different.snps[informative.families]
  family.weights  <- family.weights/sum(family.weights)

  ### compute weighted difference vectors for cases vs complements ###
  dif.vecs <- as.matrix(family.weights*cases.minus.complements[informative.families, target.snps])

  ### take the sum of the case - complement difference vectors over families ###
  sum.dif.vecs <- colSums(dif.vecs)

  ### fitness score is squared vector length of the sum of weighted difference vectors ###
  fitness.score <- sum(colSums(dif.vecs)^2)
  return(fitness.score)

}



