#' A function to assign a fitness score to a chromosome
#'
#' This function assigns a fitness score to a chromosome.
#'
#' @param case.genetic.data The genetic data of the disease affected children from case-parent trios or affected/unaffected sibling pairs. Columns are SNP allele counts, and rows are individuals.
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
#'  to be in LD if they are located on the same biological chromosome. If investigating maternal effects, where SNPs are being used as a
#'  proxy for a prenatal exposure, every entry of \code{block.ld.mat} should be set to TRUE.
#' @param weight.lookup A vector that maps a family weight to the weighted sum of the number of different SNPs and SNPs both equal to one.
#' @param case2.mat A logical matrix indicating whether, for each SNP, the case carries 2 copies of the minor allele.
#' @param case0.mat A logical matrix indicating whether, for each SNP, the case carries 0 copies of the minor allele.
#' @param n.different.snps.weight The number by which the number of different SNPs between a case and complement or unaffected sibling
#'  is multiplied in computing the family weights. Defaults to 2.
#' @param n.both.one.weight The number by which the number of SNPs equal to 1 in both the case and complement or unaffected sibling
#'  is multiplied in computing the family weights. Defaults to 1.
#' @param recode.threshold For a given SNP, the minimum test statistic required to recode and recompute the fitness score using recessive coding. Defaults to 3.
#' See the GADGETS paper for specific details.
#' @param epi.test A logical indicating whether the function should return the information required to run function \code{epistasis.test}.
#' for a given SNP. See the GADGETS paper for specific details on the implementation of this argument.
#' @return A list:
#' \describe{
#'  \item{fitness.score}{The chromosome fitness score.}
#'  \item{sum.dif.vecs}{The weighted mean difference vector corresponding to the chromosome.
#'  The magnitudes of these values are not particularly important, but the sign is useful.
#'  A positive value for a given SNP indicates the minor allele is positively associated with
#'  disease status, while a negative value implies the reference (‘wild type’) allele is
#'  positively associated with the disease.}
#'  \item{rr}{The fraction of provisional risk alleles carried by cases with the full risk set
#'  over the total number of risk alleles carried by either a case or complement with the full risk set.}
#'  \item{pseudo.t2}{The pseudo T^2 value for the chromosome.}
#'  \item{risk.set.alleles}{A vector indicating the number risk alleles a case or complement must have
#'   for each SNP in \code{target.snps} for the case or complement to be classified as having the
#'   proposed risk set. '1+' indicates at least one copy of the risk allele is required, while '2'
#'   indicates 2 copies are needed. The risk allele can be determined based on the signs of the elements
#'   of \code{sum.dif.vecs}, where a negative value indicates the major allele for a given SNP is
#'   the risk allele, while a positive value implicates the minor allele.}
#'   \item{inf.families}{An integer vector of the informative family rows. Only returned if \code{epi.test} = TRUE.}
#' }
#'
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' comp <- mom + dad - case
#' case.comp.diff <- case != comp
#' case.minus.comp <- sign(case - comp)
#' both.one.mat <- case == 1 & comp == 1
#' case2.mat <- case == 2
#' case0.mat <- case == 0
#' library(Matrix)
#' block.ld.mat <- as.matrix(bdiag(list(matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25))))
#' weight.lookup <- vapply(seq_len(6), function(x) 2^x, 1)
#' chrom.fitness.score(case, comp, case.comp.diff, c(1, 4, 7),
#'                     case.minus.comp, both.one.mat,
#'                     block.ld.mat, weight.lookup,
#'                     case2.mat, case0.mat)
#'
#' @export

chrom.fitness.score <- function(case.genetic.data, complement.genetic.data, case.comp.differences,
                                target.snps, cases.minus.complements, both.one.mat,
                                block.ld.mat, weight.lookup, case2.mat, case0.mat,
                                n.different.snps.weight = 2, n.both.one.weight = 1,
                                recode.threshold = 3, epi.test = FALSE) {

  ### pick out the differences for the target snps ###
  case.comp.diff <- case.comp.differences[, target.snps]
  case.comp.diff <- t(case.comp.diff)

  ### determine whether families are informative for the set of target.snps ###
  total.different.snps <- colSums(case.comp.diff)
  informative.families <- total.different.snps != 0
  inf.family.rows <- which(informative.families)
  n.informative.families <- sum(informative.families)
  case.inf <- t(case.genetic.data[informative.families, target.snps])
  comp.inf <- t(complement.genetic.data[informative.families, target.snps])
  cases.minus.complements.inf <- t(cases.minus.complements[informative.families, target.snps])

  ### compute weights ###
  both.one.inf <- t(both.one.mat[informative.families, target.snps])
  both.one <- colSums(both.one.inf)
  weighted.informativeness <- n.both.one.weight * both.one + n.different.snps.weight * total.different.snps[informative.families]
  family.weights <- weight.lookup[weighted.informativeness]
  invsum.family.weights <- 1/sum(family.weights)

  ### compute weighted difference vectors for cases vs complements ###
  dif.vecs <- as.matrix((family.weights *invsum.family.weights)*t(cases.minus.complements.inf))

  ### take the sum of the case - complement difference vectors over families ###
  sum.dif.vecs <- colSums(dif.vecs)

  ### determine how many cases and complements actually have the proposed risk set ###
  risk.dirs <- sign(sum.dif.vecs)
  pos.risk <- which(risk.dirs > 0)
  neg.risk <- which(risk.dirs <= 0)
  n.target <- length(target.snps)

  # Use different cases to cull out useless operations
  n.pos <- length(pos.risk)
  n.neg <- length(neg.risk)
  if (n.neg > 0) {
    n1 <- case.inf[neg.risk, , drop = FALSE] < 2
    n2 <- comp.inf[neg.risk, , drop = FALSE] < 2
  }
  if (n.pos > 0) {
    p1 <- case.inf[pos.risk, , drop = FALSE] > 0
    p2 <- comp.inf[pos.risk, , drop = FALSE] > 0
  }
  if (n.pos == 0) {
    case.high.risk <- colSums(n1) == n.target
    comp.high.risk <- colSums(n2) == n.target
  } else if (n.neg == 0) {
    case.high.risk <- colSums(p1) == n.target
    comp.high.risk <- colSums(p2) == n.target
  } else {
    case.high.risk <- (colSums(p1) + colSums(n1)) == n.target
    comp.high.risk <- (colSums(p2) + colSums(n2)) == n.target
  }
  n.case.high.risk <- sum(case.high.risk)
  n.comp.high.risk <- sum(comp.high.risk)
  case.high.inf <- case.inf[, case.high.risk, drop = FALSE]
  comp.high.inf <- comp.inf[, comp.high.risk, drop = FALSE]

  ### examine recessive mode of inheritance ####

  # initialize vector of pattern of inheritance
  risk.set.alleles <- rep("1+", length(target.snps))
  recessive <- rep(FALSE, n.target)

  if (n.case.high.risk > 0){

    # compute proportion of 2's in cases with full risk set, compare to 0.5
    if (any(pos.risk)) {

      phat <- colSums(case2.mat[which(informative.families)[case.high.risk], target.snps[pos.risk], drop = FALSE])/n.case.high.risk
      test.stat <- (phat - 0.5)/sqrt((0.5*0.5/n.case.high.risk))
      positive.recessive <- test.stat >= recode.threshold
      recessive[pos.risk] <- positive.recessive

    }

    # compute proportion of 0's in cases with full risk set, compare to 0.5
    if (any(neg.risk)) {

      phat <- colSums(case0.mat[which(informative.families)[case.high.risk], target.snps[neg.risk], drop = FALSE])/n.case.high.risk
      test.stat <- (phat - 0.5)/sqrt((0.5*0.5/n.case.high.risk))
      negative.recessive <- test.stat >= recode.threshold
      recessive[neg.risk] <- negative.recessive

    }

    ### if there are recessive SNPs, recompute the weights and associated statistics ###
    if (any(recessive)) {

      risk.set.alleles[recessive] <- "2"
      pos.outlier.cols <- pos.risk[recessive[pos.risk]]
      neg.outlier.cols <- neg.risk[recessive[neg.risk]]

      # recode instances where the model appears to be recessive
      if (length(pos.outlier.cols) > 0) {

        case.inf[pos.outlier.cols, ][case.inf[pos.outlier.cols, ] == 1] <- 0
        comp.inf[pos.outlier.cols, ][comp.inf[pos.outlier.cols, ] == 1] <- 0
        both.one.inf[pos.outlier.cols, ] <- FALSE

      }

      if (length(neg.outlier.cols) > 0) {

        case.inf[neg.outlier.cols, ][case.inf[neg.outlier.cols, ] == 1] <- 2
        comp.inf[neg.outlier.cols, ][comp.inf[neg.outlier.cols, ] == 1] <- 2
        both.one.inf[neg.outlier.cols, ] <- FALSE

      }

      # recompute the number of informative families
      cases.minus.complements <- sign(case.inf - comp.inf)
      case.comp.diff <- cases.minus.complements != 0
      total.different.snps <- colSums(case.comp.diff)
      informative.families <- total.different.snps != 0
      inf.family.rows <- inf.family.rows[informative.families]
      n.informative.families <- sum(informative.families)
      case.inf <- case.inf[, informative.families]
      comp.inf <- comp.inf[, informative.families]
      cases.minus.complements.inf <- cases.minus.complements[, informative.families]

      ### re-compute weights ###
      both.one.inf <- both.one.inf[, informative.families]
      both.one <- colSums(both.one.inf)
      weighted.informativeness <- n.both.one.weight * both.one + n.different.snps.weight * total.different.snps[informative.families]
      family.weights <- weight.lookup[weighted.informativeness]
      invsum.family.weights <- 1/sum(family.weights)

      ### compute weighted difference vectors for cases vs complements ###
      dif.vecs <- as.matrix((family.weights *invsum.family.weights)*t(cases.minus.complements.inf))
      sum.dif.vecs <- colSums(dif.vecs)

      ### re-compute proposed risk set ###
      risk.dirs <- sign(sum.dif.vecs)
      pos.risk <- which(risk.dirs > 0)
      neg.risk <- which(risk.dirs <= 0)

      n.pos <- length(pos.risk)
      n.neg <- length(neg.risk)
      if (n.neg > 0) {
        n1 <- case.inf[neg.risk, , drop = FALSE] < 2
        n2 <- comp.inf[neg.risk, , drop = FALSE] < 2
      }
      if (n.pos > 0) {
        p1 <- case.inf[pos.risk, , drop = FALSE] > 0
        p2 <- comp.inf[pos.risk, , drop = FALSE] > 0
      }
      if (n.pos == 0) {
        case.high.risk <- colSums(n1) == n.target
        comp.high.risk <- colSums(n2) == n.target
      } else if (n.neg == 0) {
        case.high.risk <- colSums(p1) == n.target
        comp.high.risk <- colSums(p2) == n.target
      } else {
        case.high.risk <- (colSums(p1) + colSums(n1)) == n.target
        comp.high.risk <- (colSums(p2) + colSums(n2)) == n.target
      }
      n.case.high.risk <- sum(case.high.risk)
      n.comp.high.risk <- sum(comp.high.risk)
      case.high.inf <- case.inf[, case.high.risk, drop = FALSE]
      comp.high.inf <- comp.inf[, comp.high.risk, drop = FALSE]

    }

  }

  ### count the number of risk alleles in those with the full risk set ###
  case.high.risk.alleles <- colSums(case.high.inf[pos.risk, , drop = FALSE]) +
    (colSums(2 - case.high.inf[neg.risk, , drop = FALSE]))
  total.case.high.risk.alleles <- sum(case.high.risk.alleles)
  comp.high.risk.alleles <- colSums(comp.high.inf[pos.risk, , drop = FALSE]) +
    (colSums(2 - comp.high.inf[neg.risk, , drop = FALSE]))
  total.comp.high.risk.alleles <- sum(comp.high.risk.alleles)

  ### compute scaling factor ###
  rr <- total.case.high.risk.alleles/(total.case.high.risk.alleles + total.comp.high.risk.alleles)
  rr <- ifelse(rr <= 0 | is.na(rr), 10^-10, rr)

  ### compute pseudo hotelling t2 ###
  mu.hat <- sum.dif.vecs
  mu.hat.mat <- matrix(rep(mu.hat, n.informative.families), nrow = n.informative.families, byrow = TRUE)
  x <- t(as.matrix(cases.minus.complements.inf))
  x.minus.mu.hat <- x - mu.hat.mat
  weighted.x.minus.mu.hat <- family.weights * x.minus.mu.hat
  cov.mat <- (invsum.family.weights) * crossprod(weighted.x.minus.mu.hat, x.minus.mu.hat)

  target.block.ld.mat <- block.ld.mat[target.snps, target.snps]
  cov.mat[!target.block.ld.mat] <- 0
  elem.vars <- sqrt(diag(cov.mat))
  sum.dif.vecs <- sum.dif.vecs/elem.vars
  zero.var <- elem.vars == 0
  sum.dif.vecs[zero.var] <- 10^-10

  # compute svd of dif.vec.cov.mat
  cov.mat.svd <- svd(cov.mat)
  cov.mat.svd$d[cov.mat.svd$d < sqrt(.Machine$double.eps)] <- 10^10

  # compute fitness score
  pseudo.t2 <- (1/(invsum.family.weights*1000)) * rowSums((t(mu.hat) %*% cov.mat.svd$u)^2/cov.mat.svd$d)
  fitness.score <- (rr^2) * pseudo.t2

  # if desired, return the required information for the epistasis test
  if (epi.test){

    high.risk.families <- inf.family.rows
    return(list(fitness.score = fitness.score, sum.dif.vecs = sum.dif.vecs, rr = rr, pseudo.t2 = pseudo.t2,
                risk.set.alleles = risk.set.alleles, inf.families = inf.family.rows))

  } else {

    return(list(fitness.score = fitness.score, sum.dif.vecs = sum.dif.vecs, rr = rr, pseudo.t2 = pseudo.t2,
                risk.set.alleles = risk.set.alleles))

  }

}

