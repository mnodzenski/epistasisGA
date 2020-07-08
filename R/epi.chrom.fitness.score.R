#' A function to assign a fitness score to a chromosome using epistasis weighting
#'
#' This function assigns a fitness score to a chromosome using epistasis weighting.
#'
#' @param case.genetic.data The genetic data of the disease affected children from case-parent trios. Columns are SNPs, and rows are individuals.
#' @param complement.genetic.data A genetic dataset from the complements of the cases, where
#' \code{complement.genetic.data} = mother SNP counts + father SNP counts - case SNP counts.
#' Columns are SNPs, rows are families.
#' @param case.comp.differences A data frame or matrix indicating \code{case.genetic.data} != \code{complement.genetic.data},
#' where rows correspond to individuals and columns correspond to snps.
#' @param target.snps A numeric vector of the columns corresponding to the collection of SNPs, or chromosome, for which the fitness score will be computed.
#' @param cases.minus.complements A matrix equal to \code{case.genetic.data} - \code{complement genetic data}.
#' @param target.snps.sum.dif.vecs A vector corresponding to the \code{sum.dif.vecs} output from \code{chrom.fitness.score} for
#' \code{target.snps}.
#' @param target.snps.risk.alleles A corresponding to the \code{risk.set.alleles} output from \code{chrom.fitness.score} for
#' \code{target.snps}.
#' @param chrom.mat A logical matrix indicating whether the SNPs in \code{case.comp.differences} are located on the same biological chromosome.
#' @param weight.lookup A vector that maps a family weight to the weighted sum of the number of different SNPs and SNPs both equal to one.
#' @param n.different.snps.weight The number by which the number of different SNPs between a case and complement is multiplied in computing the family weights. Defaults to 2.
#' @param n.both.risk.geno.weight The number by which the number of SNPs with the risk genotype in both the case and complement is multiplied in computing the family weights. Defaults to 1.
#' @param n.case.high.risk.thresh The number of cases with the provisional high risk set required to check for recessive patterns of allele inheritance.
#' @param outlier.sd The number of standard deviations from the mean allele count used to determine whether recessive allele coding is appropriate
#' @param epi.test A logical indicating whether the function should return the information required to run function \code{epi.test}.
#' for a given SNP. See the GADGET paper for specific details on the implementation of this argument.
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
#' library(Matrix)
#' chrom.mat <- as.matrix(bdiag(list(matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25))))
#' weight.lookup <- vapply(seq_len(6), function(x) 2^x, 1)
#' first.pass <- chrom.fitness.score(case, comp, case.comp.diff, c(1, 4, 7),
#'                             case.minus.comp, both.one.mat,
#'                             chrom.mat, weight.lookup)
#' epi.chrom.fitness.score(case, comp, case.comp.diff, c(1, 4, 7), case.minus.comp,
#'                         first.pass$sum.dif.vecs, first.pass$risk.set.alleles,
#'                         chrom.mat, weight.lookup)
#'
#' @export

epi.chrom.fitness.score <- function(case.genetic.data, complement.genetic.data,
                                target.snps, target.snps.sum.diff.vecs,
                                target.snps.risk.alleles, chrom.mat, weight.lookup,
                                n.different.snps.weight = 2, n.both.risk.geno.weight = 1,
                                n.case.high.risk.thresh = 20, outlier.sd = 2.5, epi.test = FALSE) {

  ### make copy of case/complement data in case originals are modified ###
  case.copy <- t(case.genetic.data[ , target.snps])
  comp.copy <- t(complement.genetic.data[ , target.snps])
  case <- case.copy
  comp <- comp.copy

  ### recode provisional recessive SNPs ###
  recessive.snps <- target.snps.risk.alleles == "2"
  if (any(recessive.snps)){

    for (rec.snp in which(recessive.snps)){

      recessive.risk.dir <- target.snps.sum.diff.vecs[rec.snp]
      if (recessive.risk.dir <= 0) {

        case[rec.snp, ][case[rec.snp, ] == 1] <- 2
        comp[rec.snp, ][comp[rec.snp, ] == 1] <- 2

      } else {

        case[rec.snp, ][case[rec.snp, ] == 1] <- 0
        comp[rec.snp, ][comp[rec.snp, ] == 1] <- 0

      }

    }

  }

  ### compute case comp differences ###
  case.comp.diff <- case != comp

  ### compute cases minus complements ###
  cases.minus.complements <- sign(case - comp)

  ### determine whether families are informative for the set of target.snps ###
  total.different.snps <- colSums(case.comp.diff)
  informative.families <- total.different.snps != 0
  inf.family.rows <- which(informative.families)
  n.informative.families <- sum(informative.families)
  case.inf <- case[, informative.families]
  comp.inf <- comp[, informative.families]
  cases.minus.complements.inf <- cases.minus.complements[, informative.families]

  ### compute weights ###
  risk.genotypes <- mapply(function(diff.vec, n.risk.alleles){

    if (sign(diff.vec) > 0){

      if (n.risk.alleles == "1+"){

        c(1,2)

      } else{

        2

      }

    } else {

      if (n.risk.alleles == "1+"){

        c(0,1)

      } else{

        0

      }

    }

  }, target.snps.sum.diff.vecs, target.snps.risk.alleles)
  both.risk.geno.mat <- matrix(FALSE, nrow = nrow(case.inf), ncol = ncol(case.inf))
  for (i in 1:ncol(both.risk.geno.mat)){

    risk.geno <- risk.genotypes[[i]]
    both.risk.geno.mat[i, ] <- case.inf[i, ] == comp.inf[i, ] & case.inf[i, ] %in% risk.geno

  }

  both.risk.geno <- colSums(both.risk.geno.mat)
  weighted.informativeness <- n.both.risk.geno.weight * both.risk.geno +
    n.different.snps.weight * total.different.snps[informative.families]
  family.weights <- weight.lookup[weighted.informativeness]
  invsum.family.weights <- 1/sum(family.weights)

  ### compute weighted difference vectors for cases vs complements ###
  dif.vecs <- as.matrix((family.weights *invsum.family.weights)*t(cases.minus.complements.inf))

  ### take the sum of the case - complement difference vectors over families ###
  sum.dif.vecs <- colSums(dif.vecs)

  ### make sure sum.dif.vecs is still consistent with proposed risk directions ###
  ### if not, change risk genotype coding to opposite allele and recompute ###
  ### sum.dif.vecs ###
  if (! all(sign(dif.vecs) == sign(target.snps.sum.diff.vecs))){

    change.these <- which(sign(dif.vecs) != sign(target.snps.sum.diff.vecs))
    for (snp in change.these){

      if (target.snps.risk.alleles[snp] == "2"){

        case[snp, ] <- case.copy[snp, ]
        comp[snp, ] <- comp.copy[snp, ]
        target.snps.risk.alleles[snp] <- "1+"

      }
      if (sign(dif.vecs)[snp] <= 0){

        new.risk.geno <- c(0,1)

      } else {

        new.risk.geno <- c(1,2)

      }
      risk.genotypes[[snp]] <- new.risk.geno

    }

    ### compute case comp differences ###
    case.comp.diff <- case != comp

    ### compute cases minus complements ###
    cases.minus.complements <- sign(case - comp)

    ### determine whether families are informative for the set of target.snps ###
    total.different.snps <- colSums(case.comp.diff)
    informative.families <- total.different.snps != 0
    inf.family.rows <- which(informative.families)
    n.informative.families <- sum(informative.families)
    case.inf <- case[informative.families, ]
    comp.inf <- comp[informative.families, ]
    cases.minus.complements.inf <- cases.minus.complements[informative.families, ]

    both.risk.geno.mat <- matrix(FALSE, nrow = nrow(case.inf), ncol = ncol(case.inf))
    for (i in 1:ncol(both.risk.geno.mat)){

      risk.geno <- risk.genotypes[[i]]
      both.risk.geno.mat[i, ] <- case.inf[i, ] == comp.inf[i, ] & case.inf[i, ] %in% risk.geno

    }

    both.risk.geno <- colSums(both.risk.geno.mat)
    weighted.informativeness <- n.both.risk.geno.weight * both.risk.geno +
      n.different.snps.weight * total.different.snps[informative.families]
    family.weights <- weight.lookup[weighted.informativeness]
    invsum.family.weights <- 1/sum(family.weights)

    ### compute weighted difference vectors for cases vs complements ###
    dif.vecs <- as.matrix((family.weights *invsum.family.weights)*t(cases.minus.complements.inf))

    ### take the sum of the case - complement difference vectors over families ###
    sum.dif.vecs <- colSums(dif.vecs)

  }

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

  ### pick out misclassifications via outlier detection, indicating recessive mode of inheritance ####

  # initialize vector of pattern of inheritance
  risk.set.alleles <- target.snps.risk.alleles
  already.recessive <- risk.set.alleles == "2"

  # only applies if we have at least 20 high risk case
  if (n.case.high.risk > n.case.high.risk.thresh) {

    case.high.risk.means <- colMeans(t(case.high.inf))
    case.high.risk.sd <- colSds(t(as.matrix(case.high.inf)))
    case.high.risk.sd[is.na(case.high.risk.sd)] <- 0
    high.outlier.thresh <- case.high.risk.means + outlier.sd * case.high.risk.sd
    low.outlier.thresh <- case.high.risk.means - outlier.sd * case.high.risk.sd
    all.high.risk <- t(cbind(case.high.inf, comp.high.inf))
    n.high.risk <- nrow(all.high.risk)
    outliers <- rep(FALSE, n.target)
    if (any(pos.risk)) {

      these.pos <- pos.risk & !already.recessive
      positive.high.risk <- all.high.risk[, these.pos, drop = FALSE]
      low.outlier.mat <- matrix(rep(low.outlier.thresh[these.pos], n.high.risk), nrow = n.high.risk,
                                byrow = TRUE)
      positive.outliers <- colSums(positive.high.risk < low.outlier.mat) > 0
      outliers[these.pos] <- positive.outliers

    }

    if (any(neg.risk)) {

      these.neg <- neg.risk & !already.recessive
      negative.high.risk <- all.high.risk[, these.neg, drop = FALSE]
      high.outlier.mat <- matrix(rep(high.outlier.thresh[these.neg], n.high.risk), nrow = n.high.risk,
                                 byrow = TRUE)
      negative.outliers <- colSums(negative.high.risk > high.outlier.mat) > 0
      outliers[these.neg] <- negative.outliers

    }

    ### if there are outliers, recompute the weights and associated statistics ###
    if (any(outliers)) {

      risk.set.alleles[outliers] <- "2"
      pos.outlier.cols <- pos.risk[outliers[pos.risk]]
      neg.outlier.cols <- neg.risk[outliers[neg.risk]]

      # recode instances where the model appears to be recessive
      if (length(pos.outlier.cols) > 0) {

        case.inf[pos.outlier.cols, ][case.inf[pos.outlier.cols, ] == 1] <- 0
        comp.inf[pos.outlier.cols, ][comp.inf[pos.outlier.cols, ] == 1] <- 0
        risk.genotypes[pos.outlier.cols] <- 2

      }

      if (length(neg.outlier.cols) > 0) {

        case.inf[neg.outlier.cols, ][case.inf[neg.outlier.cols, ] == 1] <- 2
        comp.inf[neg.outlier.cols, ][comp.inf[neg.outlier.cols, ] == 1] <- 2
        risk.genotypes[neg.outlier.cols] <- 0

      }

      ### recompute case comp differences ###
      case.comp.diff <- case != comp

      ### recompute cases minus complements ###
      cases.minus.complements <- sign(case - comp)

      ### determine whether families are informative for the set of target.snps ###
      total.different.snps <- colSums(case.comp.diff)
      informative.families <- total.different.snps != 0
      inf.family.rows <- which(informative.families)
      n.informative.families <- sum(informative.families)
      case.inf <- case[informative.families, ]
      comp.inf <- comp[informative.families, ]
      cases.minus.complements.inf <- cases.minus.complements[informative.families, ]

      both.risk.geno.mat <- matrix(FALSE, nrow = nrow(case.inf), ncol = ncol(case.inf))
      for (i in 1:ncol(both.risk.geno.mat)){

        risk.geno <- risk.genotypes[[i]]
        both.risk.geno.mat[i, ] <- case.inf[i, ] == comp.inf[i, ] & case.inf[i, ] %in% risk.geno

      }

      both.risk.geno <- colSums(both.risk.geno.mat)
      weighted.informativeness <- n.both.risk.geno.weight * both.risk.geno +
        n.different.snps.weight * total.different.snps[informative.families]
      family.weights <- weight.lookup[weighted.informativeness]
      invsum.family.weights <- 1/sum(family.weights)

      ### compute weighted difference vectors for cases vs complements ###
      dif.vecs <- as.matrix((family.weights *invsum.family.weights)*t(cases.minus.complements.inf))

      ### take the sum of the case - complement difference vectors over families ###
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

  target.chrom.mat <- chrom.mat[target.snps, target.snps]
  cov.mat[!target.chrom.mat] <- 0
  sum.dif.vecs <- sum.dif.vecs/sqrt(diag(cov.mat))

  # compute svd of dif.vec.cov.mat
  cov.mat.svd <- svd(cov.mat)
  cov.mat.svd$d[cov.mat.svd$d < sqrt(.Machine$double.eps)] <- 10^10

  # compute fitness score
  pseudo.t2 <- (1/(invsum.family.weights*1000)) * rowSums((t(mu.hat) %*% cov.mat.svd$u)^2/cov.mat.svd$d)
  fitness.score <- (rr^2) * pseudo.t2

  # if desired, return the required information for the epistasis test
  if (epi.test){

    if (any(case.high.risk) & any(comp.high.risk)){

       high.risk.families <- list(cases = inf.family.rows[case.high.risk],
                                  complements = inf.family.rows[comp.high.risk])

    } else if (any(case.high.risk) & !any(comp.high.risk)){

      high.risk.families <- list(cases = inf.family.rows[case.high.risk])

    } else if (!any(case.high.risk) & any(comp.high.risk)){

      high.risk.families <- list(complements = inf.family.rows[comp.high.risk])

    }

    return(list(fitness.score = fitness.score, sum.dif.vecs = sum.dif.vecs, rr = rr, pseudo.t2 = pseudo.t2,
                risk.set.alleles = risk.set.alleles, high.risk.families = high.risk.families))

  } else {

    return(list(fitness.score = fitness.score, sum.dif.vecs = sum.dif.vecs, rr = rr, pseudo.t2 = pseudo.t2,
                risk.set.alleles = risk.set.alleles))

  }

}

