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
#' @param min.n.risk.set A scalar indicating the minimum number of individuals whose case - control difference vector must have sign consistent with the sign of the weighted sum of the differences vectors across families. Defaults to 10.
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
                                n.different.snps.weight = 2, n.both.one.weight = 1, weight.function = identity, min.n.risk.set = 10){

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
  #family.weights  <- family.weights/sum(family.weights)
  sum.family.weights <- sum(family.weights)

  ### compute weighted difference vectors for cases vs complements ###
  dif.vecs <- as.matrix(family.weights*cases.minus.complements.inf)/sum.family.weights
  #dif.vecs <- as.matrix(family.weights*cases.minus.complements[informative.families, target.snps])

  ### take the sum of the case - complement difference vectors over families ###
  sum.dif.vecs <- colSums(dif.vecs)
  #print(sum.dif.vecs)

  ### determine how many cases and complements actually have the proposed risk set ###
  risk.dirs <- sign(sum.dif.vecs)
  pos.risk <- which(risk.dirs > 0)
  neg.risk <- which(risk.dirs <= 0)

  case <- case.genetic.data[ , target.snps]
  comp <- complement.genetic.data[ , target.snps]

  n.target <- length(target.snps)
  case.high.risk <- (rowSums(case.inf[ , pos.risk, drop = F] > 0) +  rowSums(case.inf[ , neg.risk, drop = F] < 2)) == n.target
  comp.high.risk <- (rowSums(comp.inf[ , pos.risk, drop = F] > 0) +  rowSums(comp.inf[ , neg.risk, drop = F] < 2)) == n.target

  n.case.risk <- sum(family.weights[case.high.risk])
  n.comp.risk <- sum(family.weights[comp.high.risk])
  rr <- n.case.risk/(n.case.risk  + n.comp.risk)
  rr <- ifelse(rr == 0 | is.na(rr), 10^-10, rr)
  #n.risk.ratio <- ifelse(is.na(rr), 0, rr)
  #ppv <-ifelse(is.na(rr), 0, n.case.risk/(n.case.risk + n.comp.risk))
  #n.both.risk <- sum(case.high.risk & comp.high.risk)

  #print("Case Risk Set:")
  #print(n.case.risk)
  #print("Comp Risk Set:")
  #print(n.comp.risk)
  #print("Both High Risk:")
  #print(n.both.risk)

  #risk.set.sign.mat <- matrix(rep(sign(sum.dif.vecs), n.informative.families), nrow = n.informative.families, byrow = T)
  #target.snp.signs <- sign(cases.minus.complements[informative.families, target.snps])
  #case.risk.set <- sum(rowSums(target.snp.signs == risk.set.sign.mat | target.snp.signs == 0) == ncol(risk.set.sign.mat))
  #print("N Risk Set:")
  #print(n.risk.set)
  #print(n.informative.families)
  #print("Prop Risk Set:")
  #print(n.risk.set/n.informative.families)
  #dot.prods <- family.weights*as.matrix(cases.minus.complements[informative.families, target.snps]) %*% (sum.dif.vecs*sum.family.weights)
  #mean.sum.cubed.dot.prods <- mean(dot.prods^3)
  #mean.dot.prod <- mean(dot.prods)
  #denominator <- ((1/(n.informative.families - 1))*sum((dot.prods - mean.dot.prod)^2))^(3/2)
  #sample.skew <- mean.sum.cubed.dot.prods/denominator
  #print("Skewness")
  #print(sample.skew)
  #hist(dot.prods)
  #pos.dot.prods <- ifelse(dot.prods > 0, 1, 0)
  #w.prop.pos.dot.prods <- sum(family.weights*pos.dot.prods)/sum.family.weights
  #z.score.num <- (w.prop.pos.dot.prods - 0.5)
  #z.score.denom <- w.prop.pos.dot.prods*(1 - w.prop.pos.dot.prods)*sum((1/family.weights)^2)
  #w.zscore <- z.score.num/sqrt(z.score.denom)
  #print("Zscore")
  #print(w.zscore)
  #prop.scale <- ifelse(w.prop.pos.dot.prods < 0.6, 10^-10, w.prop.pos.dot.prods)
  #prop.scale <- w.prop.pos.dot.prods
  #print(prop.scale)
  #print("Prop Positive Dot Prods:")
  #phat <- sum(as.vector(dot.prods) > 0)/n.informative.families
  #print(w.prop.pos.dot.prods)
  #zscore <- (phat - 0.5)/sqrt((1/n.informative.families)*phat*(1 - phat))
  #print("Zscore:")
  #print(zscore)

  ### Otherwise, return the squared length of the sum of the case - complement differences ###
  mu.hat <- sum.dif.vecs
  mu.hat.mat <- matrix(rep(mu.hat, n.informative.families), nrow = n.informative.families, byrow = T)
  x <- as.matrix(cases.minus.complements.inf)
  x.minus.mu.hat <- x - mu.hat.mat
  #x.minus.mu.hat <- x
  weighted.x.minus.mu.hat <- family.weights*x.minus.mu.hat

  #sum.sq.weights <- sum(family.weights^2)

  cov.mat <- (1/(sum.family.weights))*crossprod(weighted.x.minus.mu.hat, x.minus.mu.hat)
  #cov.mat <- (1/(1 - sum.sq.weights))*crossprod(weighted.x.minus.mu.hat, x.minus.mu.hat)
  #cov.mat <-  crossprod(weighted.x.minus.mu.hat, x.minus.mu.hat)
  #cov.mat <- (1/(n.informative.families))*crossprod(weighted.x.minus.mu.hat, x.minus.mu.hat)
  #cov.mat <- crossprod(weighted.x.minus.mu.hat, x.minus.mu.hat)

  target.chrom.mat <- chrom.mat[target.snps, target.snps]
  cov.mat[!target.chrom.mat] <- 0
  #print(cov.mat)
  sum.dif.vecs <- sum.dif.vecs/sqrt(diag(cov.mat))
  #sq.sum.dif.vecs <- sum.dif.vecs^2
  #element.contributions <- sq.sum.dif.vecs/sum(sq.sum.dif.vecs)

  #compute svd of dif.vec.cov.mat
  cov.mat.svd <- svd(cov.mat)
  cov.mat.svd$d[cov.mat.svd$d < sqrt(.Machine$double.eps)] <- 10^10

  #compute final fitness score using generalized inverse and hotelling
  ### If not enough indviduals with the risk set, give a very low fitness score ###
  #if (n.risk.set < min.n.risk.set){

   # fitness.score <- 10^-10

  #} else {


  #  fitness.score <-  prop.scale*sum.family.weights/1000*rowSums((t(mu.hat) %*% cov.mat.svd$u)^2/cov.mat.svd$d)
  #print("T2")
   #print(sum.family.weights/1000*rowSums((t(mu.hat) %*% cov.mat.svd$u)^2/cov.mat.svd$d))
  fitness.score <- rr*sum.family.weights/1000*rowSums((t(mu.hat) %*% cov.mat.svd$u)^2/cov.mat.svd$d)
  #fitness.score <- prop.scale*sum.family.weights/1000*rowSums((t(mu.hat) %*% cov.mat.svd$u)^2/cov.mat.svd$d)

  #}

  #if (n.comp.risk > 30){

   # fitness.score <- 10^-10

  #} else {

   # fitness.score <- rr

  #}

  #fitness.score <- (10^10)*rowSums((t(mu.hat) %*% cov.mat.svd$u)^2/cov.mat.svd$d)
  #sum.dif.vecs.sq <- sum.dif.vecs^2
  #squared.vec.length <- sum(sum.dif.vecs.sq)
  #fitness.score <- squared.vec.length
  #print("Scaled Fitness Score")
  #print((n.risk.set/n.informative.families)*fitness.score)

  return(list(fitness.score = fitness.score, sum.dif.vecs = sum.dif.vecs))

}

