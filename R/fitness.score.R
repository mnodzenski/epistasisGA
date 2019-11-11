#' A function to generate a fitness score for a subset of snps from a dataset of case genetic markers
#'
#' This function takes a subset of snps from an input genetic dataset and returns a fitness score using a k-nearest neighbors approach.
#'
#' @param case.genetic.data A genetic dataset from cases (for a dichotomous trait). Columns are snps, and rows are individuals.
#' @param target.snps A numeric vector of the columns corresponding to the snps for which the fitness score will be computed.
#' @param k A numeric scalar corresponding to the number of nearest neighbors required for computing the fitness score. See details for more information.
#' @param correct.thresh A numeric scalar between 0 and 1 indicating the minimum proportion of of cases among the nearest neighbors for a given individual for that individual to be considered correctly classified. See details for more information.
#'
#' @return A scalar, the fitness score for the given set of snps.
#'
#' @examples
#'
#' data(cases)
#' fitness.score(target.snps = c(1, 4, 7), case.genetic.data = cases)
#'
#' @importFrom stats dist
#' @export

fitness.score <- function(case.genetic.data, target.snps, k = 10, correct.thresh = 0.9){

  ### 1. Pick out the target snps from the genetic data ###
  cases <- case.genetic.data[ , target.snps]
  rownames(cases) <- rep(1, nrow(cases))

  ### 2. compute matrix of complements to the cases (i.e., the untransmitted allele counts) ###
  complements <- 2 - cases
  rownames(complements) <- rep(0, nrow(complements))

  ### 3. determine whether families are informative for the set of target.snps ###
  case.comp.diff <- cases - complements
  total.different.snps <- apply(case.comp.diff, 1, function(x) sum(x != 0))
  informative.families <- total.different.snps != 0
  n.informative.families <- sum(informative.families)

  ### 4. compute the informativeness weigth of each family (prop to # of different snps between case and complement) ###
  family.weights <- 2^total.different.snps[informative.families]

  ### 5. compute pairwise distances among all cases and complements using city block distance ###
  case.comp.data <- rbind(cases[informative.families, ], complements[informative.families, ])
  distances <- as.matrix(dist(case.comp.data, method = "manhattan"))

  ### 6. find the nearest neighbors and determine percentage of cases among neighbors ###
  ids <- rownames(cases)
  prop.neighbor.case.by.id <- sapply(1:n.informative.families, function(x){

    #remove the case itself, and the complement
    these <- c(x, x + n.informative.families)
    target.vec <- distances[-these , x]

    #grab the ids of all possible neighbors
    test.ids <- rownames(distances)[-these]

    #determine nearest neighbors by radius and compute the percentage of cases among nearest neighbors
    zeroes <- target.vec == 0
    ones <- target.vec == 1
    twos <- target.vec == 2
    if (sum(zeroes) >= k){

      neighbor.ids <- test.ids[zeroes]
      case.pct <- sum(neighbor.ids != "0")/length(neighbor.ids)

    } else if (sum(zeroes + ones) >= k){

      neighbor.ids <- test.ids[c(zeroes, ones)]
      case.pct <- sum(neighbor.ids != "0")/length(neighbor.ids)

    }  else if (sum(zeroes + ones + twos) >= k){

      neighbor.ids <- test.ids[c(zeroes, ones, twos)]
      case.pct <- sum(neighbor.ids != "0")/length(neighbor.ids)

    } else {

      case.pct <- 0

    }

    #set to zero if less than cutoff threshold
    if (case.pct < correct.thresh){

      case.pct <- 0

    }
    return(case.pct)

  })

  ### 7. compute chromosome fitness scores ###
  fitness.score <- as.numeric((family.weights %*% prop.neighbor.case.by.id)/sum(family.weights))
  return(fitness.score)

}



