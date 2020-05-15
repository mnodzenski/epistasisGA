#' A function to generate a fitness score for a subset of snps from a dataset of case genetic markers
#'
#' This function returns a fitness score for a set of snps, where the fitness score is the squared vector length of the weighted sum of difference vectors between cases and complements.
#'
#' @param case.genetic.data A genetic dataset from cases (for a dichotomous trait). Columns are snps, and rows are individuals.
#' @param complement.genetic.data A genetic dataset from the complements of the cases, where \code{complement.genetic.data} = mother snp counts + father snp counts - case snp counts. Columns are snps, rows are families. If not specified, \code{father.genetic.data} and \code{mother.genetic.data} must be specified.
#' @param case.comp.differences a data frame or matrix indicating case genetic data != complement genetic data, where rows correspond to individuals and columns correspond to snps.
#' @param target.snps A numeric vector of the columns corresponding to the snps for which the fitness score will be computed.
#' @param cases.minus.complements A matrix equal to case genetic data - complement genetic data.
#' @param both.one.mat A matrix whose elements indicate whether both the case and complement have one copy of the alternate allele, equal to (case.genetic.data == 1 & complement.genetic.data == 1).
#' @param n.different.snps.weight The number by which the number different snps between case and control is multiplied in computing the family weights. Defaults to 2.
#' @param n.both.one.weight The number by which the number of snps equal to 1 in both case and control is multiplied in computing the family weights. Defaults to 1.
#' @param weight.function A function which takes the weighted sum of the number of different snps and snps both equal to one as an argument, and returns a family weight. Defaults to the identity function.
#' @param chrom.mat A logical matrix indicating whether the snps in \code{case.comp.differences} belong to the same chromosome.
#' @param n.case.high.risk.thresh The number of cases with the provisional high risk set required to check for recessive patterns of allele inheritance.
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
#' library(Matrix)
#' chrom.mat <- as.matrix(bdiag(list(matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25))))
#' chrom.fitness.score(case, comp, case.comp.diff, c(1, 4, 7),
#'                     case.minus.comp, both.one.mat, chrom.mat)
#'
#' @export

chrom.fitness.score <- function(case.genetic.data, complement.genetic.data, case.comp.differences, target.snps,
                                cases.minus.complements, both.one.mat, chrom.mat,
                                n.different.snps.weight = 2, n.both.one.weight = 1, weight.function = function(x) 2^x,
                                n.case.high.risk.thresh = 20){

  ### pick out the differences for the target snps ###
  case.comp.diff <- case.comp.differences[ , target.snps]

  ### determine whether families are informative for the set of target.snps ###
  total.different.snps <-rowSums(case.comp.diff)
  informative.families <- total.different.snps != 0
  n.informative.families <- sum(informative.families)
  case.inf <- case.genetic.data[informative.families, target.snps]
  comp.inf <- complement.genetic.data[informative.families, target.snps]
  cases.minus.complements.inf <- cases.minus.complements[informative.families, target.snps]

  ### compute weights ###
  both.one.inf <- both.one.mat[informative.families , target.snps]
  both.one <- rowSums(both.one.inf)
  weighted.informativeness <- n.both.one.weight*both.one + n.different.snps.weight*total.different.snps[informative.families]
  family.weights <- weight.function(weighted.informativeness)
  sum.family.weights <- sum(family.weights)
  #chromosome.size <- length(target.snps)

  ### compute weighted difference vectors for cases vs complements ###
  dif.vecs <- as.matrix(family.weights*cases.minus.complements.inf)/sum.family.weights

  ### take the sum of the case - complement difference vectors over families ###
  sum.dif.vecs <- colSums(dif.vecs)

  ### determine how many cases and complements actually have the proposed risk set ###
  risk.dirs <- sign(sum.dif.vecs)
  pos.risk <- which(risk.dirs > 0)
  neg.risk <- which(risk.dirs <= 0)

  n.target <- length(target.snps)
  case.high.risk <- (rowSums(case.inf[ , pos.risk, drop = F] > 0) +  rowSums(case.inf[ , neg.risk, drop = F] < 2)) == n.target
  n.case.high.risk <- sum(case.high.risk)
  comp.high.risk <- (rowSums(comp.inf[ , pos.risk, drop = F] > 0) +  rowSums(comp.inf[ , neg.risk, drop = F] < 2)) == n.target
  #n.comp.high.risk <- sum(comp.high.risk)
  case.high.inf <- case.inf[case.high.risk, , drop = F]
  comp.high.inf <- comp.inf[comp.high.risk, , drop = F]

  ### pick out misclassifications via outlier detection, indicating recessive mode of inheritance ####

  #only applies if we have at least 20 high risk case
  if (n.case.high.risk > n.case.high.risk.thresh){

    case.high.risk.means <- colMeans(case.high.inf)
    case.high.risk.sd <- colSds(as.matrix(case.high.inf))
    case.high.risk.sd[is.na(case.high.risk.sd)] <- 0
    high.outlier.thresh <- case.high.risk.means + 2.5*case.high.risk.sd
    low.outlier.thresh <- case.high.risk.means - 2.5*case.high.risk.sd
    all.high.risk <- rbind(case.high.inf, comp.high.inf)
    n.high.risk <- nrow(all.high.risk)
    outliers <- rep(F, n.target)
    if (any(pos.risk)){

      positive.high.risk <- all.high.risk[ , pos.risk, drop = F]
      low.outlier.mat <- matrix(rep(low.outlier.thresh[pos.risk], n.high.risk), nrow = n.high.risk, byrow = T)
      positive.outliers <- colSums(positive.high.risk < low.outlier.mat) > 0
      outliers[pos.risk] <- positive.outliers

    }

    if (any(neg.risk)){

      negative.high.risk <- all.high.risk[ , neg.risk, drop = F]
      high.outlier.mat <- matrix(rep(high.outlier.thresh[neg.risk], n.high.risk), nrow = n.high.risk, byrow = T)
      negative.outliers <- colSums(negative.high.risk > high.outlier.mat) > 0
      outliers[neg.risk] <- negative.outliers

    }

    ### if there are outliers, recompute the weights and associated statistics ###
    if (any(outliers)){

      pos.outlier.cols <- pos.risk[outliers[pos.risk]]
      neg.outlier.cols <- neg.risk[outliers[neg.risk]]

      #recode instances where the model appears to be recessive
       if (length(pos.outlier.cols) > 0){

        case.inf[ , pos.outlier.cols][case.inf[ , pos.outlier.cols] == 1] <- 0
        comp.inf[ , pos.outlier.cols][comp.inf[ , pos.outlier.cols] == 1] <- 0
        both.one.inf[ , pos.outlier.cols] <- F

       }

      if (length(neg.outlier.cols) > 0){

        case.inf[ , neg.outlier.cols][case.inf[ , neg.outlier.cols] == 1] <- 2
        comp.inf[ , neg.outlier.cols][comp.inf[ , neg.outlier.cols] == 1] <- 2
        both.one.inf[ , neg.outlier.cols] <- F

      }

      #recompute the number of informative families
      cases.minus.complements <- sign(case.inf - comp.inf)
      case.comp.diff <- cases.minus.complements != 0
      total.different.snps <-rowSums(case.comp.diff)
      informative.families <- total.different.snps != 0
      n.informative.families <- sum(informative.families)
      case.inf <- case.inf[informative.families, ]
      comp.inf <- comp.inf[informative.families, ]
      cases.minus.complements.inf <- cases.minus.complements[informative.families, ]

      ### re-compute weights ###
      both.one.inf <- both.one.inf[informative.families , ]
      both.one <- rowSums(both.one.inf)
      weighted.informativeness <- n.both.one.weight*both.one + n.different.snps.weight*total.different.snps[informative.families]
      family.weights <- weight.function(weighted.informativeness)
      sum.family.weights <- sum(family.weights)

      ### re-compute weighted difference vectors for cases vs complements ###
      dif.vecs <- as.matrix(family.weights*cases.minus.complements.inf)/sum.family.weights
      sum.dif.vecs <- colSums(dif.vecs)

      ### re-compute proposed risk set ###
      risk.dirs <- sign(sum.dif.vecs)
      pos.risk <- which(risk.dirs > 0)
      neg.risk <- which(risk.dirs <= 0)

      case.high.risk <- (rowSums(case.inf[ , pos.risk, drop = F] > 0) +  rowSums(case.inf[ , neg.risk, drop = F] < 2)) == n.target
      comp.high.risk <- (rowSums(comp.inf[ , pos.risk, drop = F] > 0) +  rowSums(comp.inf[ , neg.risk, drop = F] < 2)) == n.target
      case.high.inf <- case.inf[case.high.risk, , drop = F]
      comp.high.inf <- comp.inf[comp.high.risk, , drop = F]

    }
  }

  ### count the number of risk alleles in those with the full risk set ###
  case.high.risk.alleles <- rowSums(case.high.inf[ , pos.risk, drop = F]) + (rowSums(2 - case.high.inf[ , neg.risk, drop = F]))
  total.case.high.risk.alleles <- sum(case.high.risk.alleles)
  comp.high.risk.alleles <- rowSums(comp.high.inf[ , pos.risk, drop = F]) + (rowSums(2 - comp.high.inf[ , neg.risk, drop = F]))
  total.comp.high.risk.alleles <- sum(comp.high.risk.alleles)

  ### compute scaling factor ###
  rr <- total.case.high.risk.alleles/(total.case.high.risk.alleles + total.comp.high.risk.alleles)
  rr <- ifelse(rr <= 0 | is.na(rr), 10^-10, rr)

  ### compute pseudo hotelling t2 ###
  mu.hat <- sum.dif.vecs
  mu.hat.mat <- matrix(rep(mu.hat, n.informative.families), nrow = n.informative.families, byrow = T)
  x <- as.matrix(cases.minus.complements.inf)
  x.minus.mu.hat <- x - mu.hat.mat
  weighted.x.minus.mu.hat <- family.weights*x.minus.mu.hat
  cov.mat <- (1/(sum.family.weights))*crossprod(weighted.x.minus.mu.hat, x.minus.mu.hat)

  target.chrom.mat <- chrom.mat[target.snps, target.snps]
  cov.mat[!target.chrom.mat] <- 0
  sum.dif.vecs <- sum.dif.vecs/sqrt(diag(cov.mat))

  #compute svd of dif.vec.cov.mat
  cov.mat.svd <- svd(cov.mat)
  cov.mat.svd$d[cov.mat.svd$d < sqrt(.Machine$double.eps)] <- 10^10

  #compute fitness score
  pseudo.t2 <- (sum.family.weights/1000)*rowSums((t(mu.hat) %*% cov.mat.svd$u)^2/cov.mat.svd$d)
  fitness.score <- (rr^2)*pseudo.t2

  return(list(fitness.score = fitness.score, sum.dif.vecs = sum.dif.vecs, rr = rr, pseudo.t2 = pseudo.t2))

}

