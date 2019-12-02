#' A function to generate a fitness score for a subset of snps from a dataset of case genetic markers
#'
#' This function takes a subset of snps from an input genetic dataset and returns a fitness score using a k-nearest neighbors approach.
#'
#' @param case.genetic.data A genetic dataset from cases (for a dichotomous trait). Columns are snps, and rows are individuals.
#' @param complement.genetic.data A genetic dataset representing the genetic complements to the cases (for a dichotomous trait). That is, these data correspond to the hypothetical pseudo sibling who inherited the parental alleles not transmitted to the case. Columns are snps, and rows are individuals.
#' @param case.comp.differences a data from indicating case.genetic.data != complement.genetic.data.
#' @param target.snps A numeric vector of the columns corresponding to the snps for which the fitness score will be computed.
#' @param dist.type A character string indicating the type of distance measurement. Type 'knn' performs k-nearest neighbors classifications, type 'paired' performs paired length classifications. Hotelling computes a Hotellings t2 statistic.
#' @param cases.minus.complements A matrix equal to case.genetic.data - complement.genetic.data. Required if dist.type = 'paired'.
#' @param both.one.mat A matrix whose elements indicate whether both the case and control have one copy of the alternate allele, equal to (case.genetic.data == 1 & complement.genetic.data == 1).
#' @param k A numeric scalar corresponding to the number of nearest neighbors required for computing the fitness score. See details for more information.
#' @param correct.thresh A numeric scalar between 0 and 1 indicating the minimum proportion of of cases among the nearest neighbors for a given individual for that individual to be considered correctly classified. See details for more information.
#' @return A scalar, the fitness score for the given set of snps.
#'
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' comp <- mom + dad - case
#' case.comp.diff <- case != comp
#' fitness.score(case, comp, case.comp.diff, target.snps = c(1, 4, 7), dist.type = "knn")
#'
#' @importFrom Rfast Dist
#' @export

fitness.score <- function(case.genetic.data, complement.genetic.data, case.comp.differences,
                          target.snps, dist.type, cases.minus.complements = NULL,
                          both.one.mat = NULL, k = 10, correct.thresh = 0.9){

  ###  Pick out the target snps from the case genetic data ###
  cases <- case.genetic.data[ , target.snps]

  ### Pick out target snps from the complement genetic data ###
  complements <- complement.genetic.data[ , target.snps]

  ### pick out the differences for the target snps ###
  case.comp.diff <- case.comp.differences[ , target.snps]

  ### determine whether families are informative for the set of target.snps ###
  total.different.snps <-rowSums(case.comp.diff)
  informative.families <- total.different.snps != 0
  n.informative.families <- sum(informative.families)

  ### distance computation for the knn approach ###
  if (dist.type == "knn"){

    ### compute the informativeness weigth of each family (prop to # of different snps between case and complement) ###
    family.weights <- 2^total.different.snps[informative.families]

    ### compute pairwise distances among all cases and complements using city block distance ###
    case.comp.data <- rbind(cases[informative.families, ], complements[informative.families, ])
    distances <- Dist(case.comp.data, method = "manhattan")[ , 1:n.informative.families]
    case.rows <- 1:n.informative.families
    complement.rows <- (n.informative.families + 1):(n.informative.families*2)

    ### find the nearest neighbors and determine percentage of cases among neighbors ###

    #compute the number of neighbors with within radii of 0, 1, 2
    zeroes <- distances == 0
    ones <- distances <= 1
    twos <- distances <= 2

    #now split up by radii, and compute the proportion of cases among the neighbors within the radius
    #start with radius = 0 (note: I'm subtracting 1 because the distances include the chromosome itself)
    case.zeroes <- colSums(zeroes[case.rows, ]) - 1
    total.zeroes <- colSums(zeroes) - 1
    case.comp.zeroes.ratio <- case.zeroes/total.zeroes

    #radius <= 1
    case.ones <- colSums(ones[case.rows, ]) - 1
    total.ones <- colSums(ones) - 1
    case.comp.ones.ratio <- case.ones/total.ones

    #radius <= 2
    case.twos <- colSums(twos[case.rows, ]) - 1
    total.twos <- colSums(twos) - 1
    case.comp.twos.ratio <- case.twos/total.twos

    #return a vector where the value is the proportion if there are at least k nearest neighbors and
    #the proportion of the nearest neighbors is greater than correct.thresh
    #and zero otherwise
    final.props <- rep(0, n.informative.families)
    final.props[total.twos >= k] <- case.comp.twos.ratio[total.twos >= k]
    final.props[total.ones >= k] <- case.comp.ones.ratio[total.ones >= k]
    final.props[total.zeroes >= k] <- case.comp.zeroes.ratio[total.zeroes >= k]
    final.props[final.props < correct.thresh] <- 0

    ### compute chromosome fitness scores ###
    fitness.score <- sum(family.weights*final.props)/sum(family.weights)

  } else if (dist.type == "paired"){

    #compute weights
    both.one <- rowSums(both.one.mat[informative.families , target.snps])
    family.weights <- both.one + 2*total.different.snps[informative.families]

    #difference vectors between cases and complements
    dif.vecs <- cases.minus.complements[informative.families, target.snps]

    #sum the difference vectors
    sum.dif.vecs <- colSums(dif.vecs)

    #also compute the average difference vector
    ave.dif.vec <- sum.dif.vecs/n.informative.families

    #determine whether the dot product between each family's difference vector
    #and the average difference vector is positive
    dot.prods <- dif.vecs %*% ave.dif.vec
    non.pos.dot.prods <- dot.prods <= 0

    #weight the vector sum, giving weight zero to families with negative dot products
    #and otherwise the family weights specified above
    family.weights[non.pos.dot.prods] <- 0
    weighted.dif.vecs <- family.weights*dif.vecs

    #squared vector length of the sum of difference vectors
    sq.length.sum.dif.vecs <- sum(colSums(weighted.dif.vecs)^2)

    #sd of the elements of the weighted difference vectors
    weighted.dif.vec.sd <- sd(colSums(weighted.dif.vecs)^2)

    #expected squared length of the sum of the difference vectors under the null
    expected.sq.length.sum.dif.vecs <- sum(weighted.dif.vecs^2)

    #fitness score as the ratio of the observed to expected squared vector length, inverse weigthed by
    #1 plus the variance of the absolute value of elements of the vector
    fitness.score <- (1/(1 + weighted.dif.vec.sd))*(sq.length.sum.dif.vecs/(expected.sq.length.sum.dif.vecs))

  } else if (dist.type == "Hotelling"){

    #difference vectors for informative families
    dif.vecs <- as.matrix(cases.minus.complements[informative.families , target.snps])

    #average difference vector
    ave.dif.vec <- colMeans(dif.vecs)

    #observed variance covariance matrix times n
    dif.vec.cov.mat <- crossprod(dif.vecs)

    #compute svd of dif.vec.cov.mat
    cov.mat.svd <- svd(dif.vec.cov.mat)

    #compute the hotelling t2 stat
    t2 <- length(dif.vecs) * rowSums((t(ave.dif.vec) %*% cov.mat.svd$u)^2/cov.mat.svd$d)
    fitness.score <- t2

  }

  return(fitness.score)

}



