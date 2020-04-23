#' A function to code the genetic data according to specific models of inheritance identify using a data driven approach.
#'
#' This function performs model selection on each input snp and recodes the data based on the results.
#'
#' @param case.genetic.data A genetic dataset from cases (for a dichotomous trait). Columns are snps, and rows are individuals.
#' @param complement.genetic.data A genetic dataset from the complements of the cases, where \code{complement.genetic.data} = mother snp counts + father snp counts - case snp counts. Columns are snps, rows are families. If not specified, \code{father.genetic.data} and \code{mother.genetic.data} must be specified.
#' @param father.genetic.data The genetic data for the father of the case. Columns are snps, rows are individuals. Does not need to be specified if \code{complement.genetic.data} is specified.
#' @param mother.genetic.data The genetic data for the mother of the case. Columns are snps, rows are individuals. Does not need to be specified if \code{complement.genetic.data} is specified.
#' @param chrom.mat A logical matrix indicating whether the snps in the input genetic data belong to the same chromosome.
#' @param min.allele.freq The minimum minor allele frequency in the parents required for a snp to be considered as a potential GA solution. Any snps with MAF < \code{min.allele.freq} in the parents will be omitted. Defaults to 0.01.
#' @param bp.param The BPPARAM argument to be passed to bplapply when running conditional logistic regressions on each input SNP. If using a cluster computer, this parameter needs to be set with care.
#' @return A list, whose first element is a matrix of case genetic data, second element is complment genetic data, and third element is a vector of chi-square statistics, fourth element is a list of column numbers in the original dataset to whom the .
#'
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' library(Matrix)
#' chrom.mat <- as.matrix(bdiag(list(matrix(rep(TRUE, 2500^2), nrow = 2500),
#'                               matrix(rep(TRUE, 2500^2), nrow = 2500),
#'                               matrix(rep(TRUE, 2500^2), nrow = 2500),
#'                               matrix(rep(TRUE, 2500^2), nrow = 2500))))
#' res <- preprocess.genetic.data(case[, 1:3], father.genetic.data = dad[ , 1:3], mother.genetic.data = mom[ , 1:3],
#'                  chrom.mat = chrom.mat[ , 1:3])
#'
#' @importFrom matrixStats colSds rowMaxs
#' @importFrom data.table data.table rbindlist setorder
#' @importFrom stats rbinom sd
#' @importFrom survival clogit
#' @importFrom BiocParallel bplapply
#' @importFrom survival clogit strata coxph Surv
#' @export

preprocess.genetic.data <- function(case.genetic.data, complement.genetic.data = NULL, father.genetic.data = NULL,
                                    mother.genetic.data = NULL, chrom.mat, min.allele.freq = 0.025,
                                    bp.param = bpparam()){

  #make sure the appropriate genetic data is included
  if (is.null(complement.genetic.data) & is.null(father.genetic.data) & is.null(mother.genetic.data)){

    stop("Must include complement.genetic.data or both father.genetic.data and mother.genetic.data")

  }

  ### find the snps with MAF < minimum threshold in the cases ###
  if (!is.null(father.genetic.data) & !is.null(mother.genetic.data)){

    alt.allele.freqs <- colSums(father.genetic.data + mother.genetic.data)/(4*nrow(father.genetic.data))
    minor.alleles <- alt.allele.freqs < 0.5
    below.maf.threshold <- alt.allele.freqs > (1 - min.allele.freq) | alt.allele.freqs < min.allele.freq
    original.col.numbers <- which(!below.maf.threshold)
    names(original.col.numbers) <- NULL

    ### Compute the complement data ###
    complement.genetic.data <- father.genetic.data + mother.genetic.data - case.genetic.data

    ###recode the case and complement data so that 1 indicates a copy of the minor allele ###
    case.genetic.data[ , !minor.alleles] <- 2 - case.genetic.data[ , !minor.alleles]
    complement.genetic.data[ , !minor.alleles] <- 2 - complement.genetic.data[ , !minor.alleles]

    ### remove the snps not meeting the required allele frequency threshold ###
    father.genetic.data <- father.genetic.data[ , !below.maf.threshold]
    mother.genetic.data <- mother.genetic.data[ , !below.maf.threshold]
    case.genetic.data <- case.genetic.data[ , !below.maf.threshold]
    chrom.mat <- chrom.mat[!below.maf.threshold , !below.maf.threshold]

  } else if (!is.null(complement.genetic.data)){

    alt.allele.freqs <- colSums(case.genetic.data + complement.genetic.data)/(4*nrow(case.genetic.data))
    minor.alleles <- alt.allele.freqs < 0.5
    below.maf.threshold <- alt.allele.freqs > (1 - min.allele.freq) | alt.allele.freqs < min.allele.freq
    original.col.numbers <- which(!below.maf.threshold)
    names(original.col.numbers) <- NULL

    ###recode the case and complement data so that 1 indicates a copy of the minor allele ###
    case.genetic.data[ , !minor.alleles] <- 2 - case.genetic.data[ , !minor.alleles]
    complement.genetic.data[ , !minor.alleles] <- 2 - complement.genetic.data[ , !minor.alleles]

    ### remove the snps not meeting the required allele frequency threshold ###
    case.genetic.data <- case.genetic.data[ , !below.maf.threshold]
    complement.genetic.data <- complement.genetic.data[ , !below.maf.threshold]
    chrom.mat <- chrom.mat[!below.maf.threshold , !below.maf.threshold]

  }

  ### use conditional logistic regression to estimate model of inherritance ###
  case.status <- c(rep(1, nrow(case.genetic.data)), rep(0, nrow(complement.genetic.data)))
  ids <- rep(1:nrow(case.genetic.data), 2)
  n <- nrow(case.genetic.data)
  res.list <- bplapply(1:ncol(case.genetic.data), function(snp){

    case.snp <- factor(case.genetic.data[ , snp], levels = c(0, 1, 2))
    comp.snp <- factor(complement.genetic.data[ , snp], levels = c(0,1,2))

    #get p-value of association from conditional logistic regression
    case.comp.geno <- as.numeric(as.character(c(case.snp, comp.snp)))
    clogit.res <- clogit(case.status ~ case.comp.geno + strata(ids), method = "approximate")
    clogit.chisq <- summary(clogit.res)$logtest[1]

    return(list(case.snp = case.snp, comp.snp = comp.snp, chisq = clogit.chisq))

  }, BPPARAM = bp.param)

  #combine results into new dataframes
  case.genetic.data <- do.call("cbind", lapply(res.list, function(x) as.numeric(as.character(x$case.snp))))
  complement.genetic.data <- do.call("cbind", lapply(res.list, function(x) as.numeric(as.character(x$comp.snp))))
  chisq.stats <- do.call("c", lapply(res.list, function(x) x$chisq))

  return(list(case.genetic.data = case.genetic.data, complement.genetic.data = complement.genetic.data, chisq.stats = chisq.stats, original.col.numbers = original.col.numbers,
              chrom.mat = chrom.mat, minor.allele.vec = minor.alleles))


}
