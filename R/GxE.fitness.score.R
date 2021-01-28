#' A function to assign a fitness score to a chromosome
#'
#' This function assigns a fitness score to a chromosome.
#'
#' @param case.genetic.data The genetic data of the disease affected children from case-parent trios or affected/unaffected sibling pairs.
#'  Columns are SNP allele counts, and rows are individuals.
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
#'  to be in LD if they are located on the same biological chromosome.
#' @param weight.lookup A vector that maps a family weight to the weighted sum of the number of different SNPs and SNPs both equal to one.
#' @param case2.mat A logical matrix indicating whether, for each SNP, the case carries 2 copies of the minor allele.
#' @param case0.mat A logical matrix indicating whether, for each SNP, the case carries 0 copies of the minor allele.
#' @param n.different.snps.weight The number by which the number of different SNPs between a case and complement/unaffected sibling
#'  is multiplied in computing the family weights. Defaults to 2.
#' @param n.both.one.weight The number by which the number of SNPs equal to 1 in both the case and complement/unaffected sibling
#'  is multiplied in computing the family weights. Defaults to 1.
#' @param recode.threshold For a given SNP, the minimum test statistic required to recode and recompute the fitness score using recessive coding. Defaults to 3.
#' See the GADGETS paper for specific details.
#' @param exposure A categorical factor vector corresponding to environmental exposures of the cases.
#' @return A list:
#' \describe{
#'  \item{fitness.score}{The chromosome fitness score.}
#'  \item{sum.dif.vecs}{The weighted mean difference vector corresponding to the chromosome,
#'  with each element divided by it's pseudo-standard error.
#'  The magnitudes of these values are not particularly important, but the sign is useful.
#'  A positive value for a given SNP indicates the minor allele is positively associated with
#'  disease status, while a negative value implies the reference (‘wild type’) allele is
#'  positively associated with the disease.}
#'  \item{q}{The fraction of provisional risk alleles carried by cases with the full risk set
#'  over the total number of risk alleles carried by either a case or complement with the full risk set.}
#'  \item{pseudo.t2}{The pseudo T^2 value for the chromosome.}
#'  \item{risk.set.alleles}{A vector indicating the number risk alleles a case or complement must have
#'   for each SNP in \code{target.snps} for the case or complement to be classified as having the
#'   proposed risk set. '1+' indicates at least one copy of the risk allele is required, while '2'
#'   indicates 2 copies are needed. The risk allele can be determined based on the signs of the elements
#'   of \code{sum.dif.vecs}, where a negative value indicates the major allele for a given SNP is
#'   the risk allele, while a positive value implicates the minor allele.}
#'   \item{inf.families}{An integer vector of the informative family rows. Only returned if \code{epi.test} = TRUE.}
#'   \item{xbar}{The weighted mean difference sign vector. Only returned if \code{GxE} is TRUE.}
#'   \item{sigma.hat}{The pseudo covariance matrix of the weighted mean difference sign vector. Only returned if \code{GxE} is TRUE.}
#'   \item{w}{The sum of the family weights. Only returned if \code{GxE} is TRUE.}
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
#' set.seed(11)
#' exposure <- factor(rbinom(nrow(case), 1, 0.3))
#' GxE.fitness.score(case, comp, case.comp.diff, c(1, 4, 7),
#'                     case.minus.comp, both.one.mat,
#'                     block.ld.mat, weight.lookup,
#'                     case2.mat, case0.mat, exposure)
#'
#' @export

GxE.fitness.score <- function(case.genetic.data, complement.genetic.data, case.comp.differences,
                                target.snps, cases.minus.complements, both.one.mat,
                                block.ld.mat, weight.lookup, case2.mat, case0.mat, exposure,
                                n.different.snps.weight = 2, n.both.one.weight = 1,
                                recode.threshold = 3) {

  ### divide the input data based on exposure and get components for fitness score ###
  score.by.exposure <- lapply(unique(exposure), function(exp.level){

    these.rows <- exposure == exp.level
    case.genetic.data <- case.genetic.data[these.rows, ]
    complement.genetic.data <- complement.genetic.data[these.rows, ]
    case.comp.differences <- case.comp.differences[these.rows, ]
    cases.minus.complements <- cases.minus.complements[these.rows, ]
    both.one.mat <- both.one.mat[these.rows, ]
    case2.mat <- case2.mat[these.rows, ]
    case0.mat <- case0.mat[these.rows, ]
    chrom.fitness.score(case.genetic.data, complement.genetic.data, case.comp.differences,
                        target.snps, cases.minus.complements, both.one.mat,
                        block.ld.mat, weight.lookup, case2.mat, case0.mat,
                        n.different.snps.weight, n.both.one.weight,
                        recode.threshold, epi.test = FALSE, GxE = TRUE)

  })

  ### compute two sample hotelling for each pairwise comparison ###
  all.pairs <- combn(length(unique(exposure)), 2)
  pair.scores.list <- lapply(seq(1, ncol(all.pairs)), function(pair.number){

    # pick out the required pieces
    these.two <- all.pairs[ , pair.number]
    exp1 <- these.two[[1]]
    exp2 <- these.two[[2]]
    exp1.list <- score.by.exposure[[exp1]]
    exp2.list <- score.by.exposure[[exp2]]

    # mean difference vector
    xbar <- exp1.list$xbar
    ybar <- exp2.list$xbar
    xbar.minus.ybar <- xbar - ybar

    # cov mat
    sigma.x <- exp1.list$sigma*(exp1.list$w)
    sigma.y <- exp2.list$sigma*(exp2.list$w)
    w1 <- exp1.list$w
    w2 <- exp2.list$w
    q1 <- exp1.list$q
    q2 <- exp2.list$q
    sigma.hat <- (sigma.x + sigma.y)/(w1 + w2)

    # svd of cov mat
    sigma.hat.svd <- svd(sigma.hat)
    sigma.hat.svd$d[sigma.hat.svd$d < sqrt(.Machine$double.eps)] <- 10^10

    # two sample hotelling stat, using adjusted mean difference
    weight.scalar <- (w1*w2)/(w1 + w2)
    adj.xbar.minus.ybar <- q1*xbar - q2*ybar
    s <- (weight.scalar/1000)*rowSums((t(adj.xbar.minus.ybar) %*% sigma.hat.svd$u)^2/sigma.hat.svd$d)

    # return two-sample hotelling, difference vectors/se's
    se <- sqrt(diag(sigma.hat))
    std.diff.vecs <- adj.xbar.minus.ybar/se
    std.diff.vecs[se == 0] <- 10^-10

    # also use the allele coding for the higher scoring set
    if (exp2.list$fitness.score >= exp1.list$fitness.score){

      risk.set.alleles <- exp2.list$risk.set.alleles
      high.risk.exposure <- unique(exposure)[exp2]
      low.risk.exposure <- unique(exposure)[exp1]
      diff.vec.signs <- sign(ybar)

    } else {

      risk.set.alleles <- exp1.list$risk.set.alleles
      high.risk.exposure <- unique(exposure)[exp1]
      low.risk.exposure <- unique(exposure)[exp2]
      diff.vec.signs <- sign(xbar)

    }

    return(list(s = s, std.diff.vecs = std.diff.vecs, diff.vec.signs = diff.vec.signs,
                risk.set.alleles = risk.set.alleles, high.risk.exposure = high.risk.exposure,
                low.risk.exposure = low.risk.exposure))

  })

  # pick out the biggest difference
  pair.scores <- vapply(pair.scores.list, function(x) x$s, 1.0)
  largest.stat <- which(pair.scores == max(pair.scores))

  # return the largest pairwise stat as the fitness score
  fitness.score <- pair.scores[largest.stat]
  std.diff.vecs <- pair.scores.list[[largest.stat]]$std.diff.vecs
  diff.vec.signs <- pair.scores.list[[largest.stat]]$diff.vec.signs
  risk.set.alleles <- pair.scores.list[[largest.stat]]$risk.set.alleles
  high.risk.exposure <- pair.scores.list[[largest.stat]]$high.risk.exposure
  low.risk.exposure <- pair.scores.list[[largest.stat]]$low.risk.exposure


  return(list(fitness.score = fitness.score, sum.dif.vecs = std.diff.vecs, risk.set.signs = diff.vec.signs,
              risk.set.alleles = risk.set.alleles, high.risk.exposure = high.risk.exposure,
              low.risk.exposure = low.risk.exposure))

}

