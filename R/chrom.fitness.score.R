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
                                n.different.snps.weight = 2, n.both.one.weight = 1, weight.function = identity, min.n.risk.set = 10){

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

  ### determine how many cases actually have the proposed risk set ###
  risk.set.sign.mat <- matrix(rep(sign(sum.dif.vecs), n.informative.families), nrow = n.informative.families, byrow = T)
  target.snp.signs <- sign(case.comp.diff[informative.families, ])
  n.risk.set <- sum(rowSums(risk.set.sign.mat == target.snp.signs) == ncol(risk.set.sign.mat))

  ### If not enough indviduals with the risk set, give a very low fitness score ###
  if (n.risk.set < min.n.risk.set){

    fitness.score <- 10^-10

  } else {

  ### Otherwise, return the squared length of the sum of the case - complement differences ###
    fitness.score <- sum(sum.dif.vecs^2)

  }

  ### compute the average difference vector ###
  #ave.dif.vec <- sum.dif.vecs/n.informative.families

  ### compute dot product of difference vectors with the average difference vector ###
  #dot.prods <- as.numeric(dif.vecs %*% ave.dif.vec)
  #ave.dif.vec.mat <- matrix(rep(sign(ave.dif.vec), nrow(dif.vecs)), byrow = T, nrow = nrow(dif.vecs))
  #sign.comp.mat <- sign(dif.vecs) == ave.dif.vec.mat
  #keep.these <- rowSums(sign.comp.mat) == length(ave.dif.vec)

  ### keep the families with a positive dot product ###
  #keep.these <- dot.prods > 0

  ### get final difference vectors ###
  #if( any(keep.these)){

 #   sum.dif.vecs <- colSums(dif.vecs[keep.these, , drop = F])

    ### fitness score is squared vector length of the sum of weighted difference vectors ###
#    fitness.score <- sum(sum.dif.vecs^2)

#  } else{

 #   sum.dif.vecs <- rep(10^-6, length(target.snps))
#    fitness.score <- 0

#  }

  return(list(fitness.score = fitness.score, sum.dif.vecs = sum.dif.vecs))

}



