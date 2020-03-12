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
#' chrom.mat <- as.matrix(bdiag(list(matrix(rep(TRUE, 2500^2), nrow = 2500),
#'                               matrix(rep(TRUE, 2500^2), nrow = 2500),
#'                               matrix(rep(TRUE, 2500^2), nrow = 2500),
#'                               matrix(rep(TRUE, 2500^2), nrow = 2500))))
#' chrom.fitness.score(case.comp.diff, target.snps = c(1, 4, 7), case.minus.comp, both.one.mat, chrom.mat)
#'
#' @export

chrom.fitness.score <- function(case.genetic.data, complement.genetic.data, case.comp.differences, target.snps, cases.minus.complements, both.one.mat, chrom.mat,
                                n.different.snps.weight = 2, n.both.one.weight = 1, weight.function = identity){

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
  both.one <- rowSums(both.one.mat[informative.families , target.snps])
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
  n.comp.high.risk <- sum(comp.high.risk)
  case.high.inf <- case.inf[case.high.risk, , drop = F]
  comp.high.inf <- comp.inf[comp.high.risk, , drop = F]

  ### pick out misclassifications via outlier detection ####
  if (n.case.high.risk > 0){

    case.high.risk.means <- colMeans(case.high.inf)
    case.high.risk.sd <- colSds(as.matrix(case.high.inf))
    high.outlier.thresh <- case.high.risk.means + 2.5*case.high.risk.sd
    low.outlier.thresh <- case.high.risk.means - 2.5*case.high.risk.sd
    if(any(neg.risk)){

      neg.risk.thresh <- high.outlier.thresh[neg.risk]
      case.neg.test <- case.high.inf[ , neg.risk, drop = F]
      if (n.comp.high.risk > 0){

        comp.neg.test <- comp.high.inf[ , neg.risk, drop = F]

      }

    } else {

      neg.risk.thresh <- 1
      case.neg.test <- matrix(rep(2, n.case.high.risk), ncol= 1)
      comp.neg.test <- matrix(rep(2, n.comp.high.risk), ncol= 1)

    }
    if (any(pos.risk)){

      pos.risk.thresh <- low.outlier.thresh[pos.risk]
      case.pos.test <- case.high.inf[, pos.risk, drop = F]
      if (n.comp.high.risk > 0){

        comp.pos.test <- comp.high.inf[, pos.risk, drop = F]

      }

    } else {

      pos.risk.thresh <- 1
      case.pos.test <- matrix(rep(0, n.case.high.risk), ncol= 1)
      comp.pos.test <- matrix(rep(0, n.comp.high.risk), ncol= 1)

    }
    n.rows <- max(n.case.high.risk, n.comp.high.risk)
    high.outlier.mat <- matrix(rep(neg.risk.thresh, n.rows), nrow = n.rows, byrow = T)
    low.outlier.mat <- matrix(rep(pos.risk.thresh, n.rows), nrow = n.rows, byrow = T)

    case.outliers <- (rowSums(case.pos.test >= low.outlier.mat[1:n.case.high.risk, , drop = F]) +
                        rowSums(case.neg.test <= high.outlier.mat[1:n.case.high.risk, , drop = F])) < n.target
    case.high.inf <- case.high.inf[!case.outliers, , drop = F]

    ### count the number of risk alleles in those with the full risk set ###
    case.high.risk.alleles <- rowSums(case.high.inf[ , pos.risk, drop = F]) + (rowSums(2 - case.high.inf[ , neg.risk, drop = F]))
    total.case.high.risk.alleles <- sum(case.high.risk.alleles)

    if (n.comp.high.risk > 0){

      comp.outliers <- (rowSums(comp.pos.test >= low.outlier.mat[1:n.comp.high.risk, , drop = F]) +
                          rowSums(comp.neg.test <= high.outlier.mat[1:n.comp.high.risk, , drop = F])) < n.target
      comp.high.inf <- comp.high.inf[!comp.outliers, , drop = F]

      ### count the number of risk alleles in those with the full risk set ###
      comp.high.risk.alleles <- rowSums(comp.high.inf[ , pos.risk, drop = F]) + (rowSums(2 - comp.high.inf[ , neg.risk, drop = F]))
      total.comp.high.risk.alleles <- sum(comp.high.risk.alleles)

      ### compute scaling factor ###
      rr <- total.case.high.risk.alleles/(total.case.high.risk.alleles + total.comp.high.risk.alleles)

    } else {

      rr <- 1

    }

    ### compute scaling factor ###
    #n.case.risk <- sum(family.weights[case.high.risk])
    #n.comp.risk <- sum(family.weights[comp.high.risk])
    #rr <- n.case.risk/(n.case.risk  + n.comp.risk)
    #rr <- sum(pooled.dif.vec)/sqrt(sum(pooled.dif.vec^2))
    #rr <- ifelse(is.na(rr), 10^-10, weight.function(rr))
    #rr <- sum(case.high.risk)/(sum(case.high.risk) + sum(comp.high.risk))
    #rr <- sum(case.high.risk)/(sum(case.high.risk) + sum(comp.high.risk))
    #rr <- sum(case.high.risk.alleles)/(sum(case.high.risk.alleles) + sum(comp.high.risk.alleles))

  } else {

    rr <- 10^-10

  }

  rr <- ifelse(rr <= 0 | is.na(rr), 10^-10, rr)
  #print(rr)

  ### Otherwise, return the squared length of the sum of the case - complement differences ###
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
  #fitness.score <- (rr)*(sum.family.weights/1000)*rowSums((t(mu.hat) %*% cov.mat.svd$u)^2/cov.mat.svd$d)
  pseudo.t2 <- (sum.family.weights/1000)*rowSums((t(mu.hat) %*% cov.mat.svd$u)^2/cov.mat.svd$d)
  fitness.score <- (rr^2)*pseudo.t2
  #fitness.score <- (rr)* sqrt(((sum.family.weights - chromosome.size)/(chromosome.size*(sum.family.weights - 1)))*sum.family.weights*rowSums((t(mu.hat) %*% cov.mat.svd$u)^2/cov.mat.svd$d))

  return(list(fitness.score = fitness.score, sum.dif.vecs = sum.dif.vecs, rr = rr, pseudo.t2 = pseudo.t2))

}

