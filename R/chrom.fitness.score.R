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
#' @param weight.function A function which takes the weighted sum of the number of different snps and snps both equal to one as an argument, and returns a family weight. Defaults to the identity function.
#' @param min.n.risk.set A scalar indicating the minimum number of individuals whose case - control difference vector must have sign consistent with the sign of the weighted sum of the differences vectors across families. Defaults to 10.
#' @param ld.mat A matrix of ld estimates among snps for the cases. Defaults to NULL.
#' @param max.ld A numeric indicating the maximimum ld value in the \code{ld.mat} matrix allowed among elements of a chromosome. If any elements of the chromosome have ld estimates above this threshold, the score will be set to 10^-10. Defaults to NULL, meaning there are no ld restrictions.
#' @return A list whose first element is the fitness score and second element is the sum of weighted difference vectors for the target snps.
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
                                n.different.snps.weight = 2, n.both.one.weight = 1, weight.function = identity, min.n.risk.set = 10,
                                ld.mat = NULL, max.ld = NULL){

  ### pick out the differences for the target snps ###
  case.comp.diff <- case.comp.differences[ , target.snps]

  ### determine whether families are informative for the set of target.snps ###
  total.different.snps <-rowSums(case.comp.diff)
  informative.families <- total.different.snps != 0
  n.informative.families <- sum(informative.families)

  ### compute weights ###
  both.one <- rowSums(both.one.mat[informative.families , target.snps])
  weighted.informativeness <- n.both.one.weight*both.one + n.different.snps.weight*total.different.snps[informative.families]
  family.weights <- weight.function(weighted.informativeness)
  family.weights  <- family.weights/sum(family.weights)

  ### compute weighted difference vectors for cases vs complements ###
  dif.vecs <- as.matrix(family.weights*cases.minus.complements[informative.families, target.snps])

  ### take the sum of the case - complement difference vectors over families ###
  sum.dif.vecs <- colSums(dif.vecs)

  ### if there are ld restrictions, determine whether the score should be set to 10^-10 ###
  if (!is.null(max.ld)){

    target.ld.mat <- ld.mat[ target.snps, target.snps]
    ld.vec <- target.ld.mat[upper.tri(target.ld.mat)]
    max.obs.ld <- max(ld.vec)

  #otherwise just give values that allow the algorithm to proceed to the score
  } else {

    max.obs.ld <- 0
    max.ld <- 1

  }

  if (max.obs.ld >= max.ld){

      fitness.score <- 10^-10

  } else {

    ### determine how many cases actually have the proposed risk set ###
    risk.set.sign.mat <- matrix(rep(sign(sum.dif.vecs), n.informative.families), nrow = n.informative.families, byrow = T)
    target.snp.signs <- sign(case.comp.diff[informative.families, ])
    n.risk.set <- sum(rowSums(risk.set.sign.mat == target.snp.signs) == ncol(risk.set.sign.mat))

    ### If not enough indviduals with the risk set, give a very low fitness score ###
    if (n.risk.set < min.n.risk.set){

      fitness.score <- 10^-10

    } else {

      ### Otherwise, return the squared length of the sum of the case - complement differences ###
      mean.diff.vec <- (1/n.informative.families)*sum.dif.vecs
      mean.diff.mat <- matrix(rep(mean.diff.vec, n.informative.families), nrow = n.informative.families, byrow = T)
      cov.mat <- (1/(n.informative.families - 1))*crossprod(dif.vecs - mean.diff.mat)

      #compute svd of dif.vec.cov.mat
      cov.mat.svd <- svd(cov.mat)

      #compute final fitness score using generalized inverse and hotelling
      fitness.score <- n.informative.families*rowSums((t(mean.dif.vec) %*% cov.mat.svd$u)^2/cov.mat.svd$d)
      #fitness.score <- sum(sum.dif.vecs^2)

    }

  }
  return(list(fitness.score = fitness.score, sum.dif.vecs = sum.dif.vecs))

}



