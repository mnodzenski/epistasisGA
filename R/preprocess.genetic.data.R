#' A function to pre-process case-parent triad of affected/unaffected sibling data.
#'
#' This function performs several pre-processing steps, intended for use before function run.gadgets.
#'
#' @param case.genetic.data The genetic data of the disease affected children from case-parent trios or affected/unaffected sibling pairs. Columns are SNP allele counts, and rows are individuals.
#' The ordering of the columns must be consistent with the LD structure specified in \code{block.ld.mat}.
#' @param complement.genetic.data A genetic dataset from the complements of the cases, where
#' \code{complement.genetic.data} = mother SNP counts + father SNP counts - case SNP counts. If using affected/unaffected siblings
#' this should be the genotypes for the unaffected siblings.
#' Columns are SNP allele counts, rows are families. If not specified, \code{father.genetic.data} and \code{mother.genetic.data} must be specified.
#' @param father.genetic.data The genetic data for the fathers of the cases. Columns are SNP allele counts, rows are individuals.
#' Does not need to be specified if \code{complement.genetic.data} is specified.
#' @param mother.genetic.data The genetic data for the mothers of the cases. Columns are SNP allele counts, rows are individuals.
#' Does not need to be specified if \code{complement.genetic.data} is specified.
#' @param block.ld.mat A logical, block diagonal matrix indicating whether the SNPs in \code{case.genetic.data} should be considered
#'  to be in linkage disequilibrium. Note that this means the ordering of the columns (SNPs) in \code{case.genetic.data} must be consistent
#'  with the LD blocks specified in \code{ld.block.mat}. In the absence of outside information, a reasonable default is to consider SNPs
#'  to be in LD if they are located on the same biological chromosome.
#' @param min.allele.freq The minimum minor allele frequency required for a SNP to be considered for inclusion in the genetic algorithm.
#' Any SNPs with MAF < \code{min.allele.freq} in the parents, or the combined group of affected and unaffected siblings, will be filtered out. Defaults to 0 (no filtering).
#' @param bp.param The BPPARAM argument to be passed to bplapply when estimating marginal disease associations for each SNP.
#'  If using a cluster computer, this parameter needs to be set with care. See \code{BiocParallel::bplapply} for more details
#' @param snp.sampling.probs A vector indicating the sampling probabilities of the SNPs in \code{case.genetic.data}. SNPs will be sampled in the
#' genetic algorithm proportional to the values specified. If not specified, by default, chi-square statistics of association will be computed for
#' each SNP, and sampling will be proportional to the root of these statistics. If user specified, the  vector values need not sum to 1, they just need to be positive
#' real numbers. See argument \code{prob} from function \code{sample} for more details.
#' @return A list containing the following:
#' \describe{
#'  \item{case.genetic.data}{The pre-processed version of the case genetic data. Any missing genotypes for a given family will be coded as -9 for both case and complement,
#'  resulting in that family being uninformative for that SNP.}
#'  \item{complement.genetic.data}{Pre-processed complement or unaffected sibling genetic data. If mother and father data are input,
#'  the complement genetic data are first created and then pre-processed. Any missing genotypes for a given family will be coded as -9 for both case and complement,
#'  resulting in that family being uninformative for that SNP.}
#'  \item{chisq.stats}{A vector of chi-square statistics corresponding to marginal SNP-disease associations, if \code{snp.sampling.probs}
#'  is not specified, and \code{snp.sampling.probs} if specified.}
#'  \item{original.col.numbers}{A vector indicating the original column number of each non-filtered SNP remaining in the analysis data.}
#'  \item{block.ld.mat}{The pre-processed version of \code{block.ld.mat}.}
#'  \item{minor.allele.vec}{A vector indicating whether the alternate allele was the minor allele for each column in the input data.}
#' }
#'
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' library(Matrix)
#' block.ld.mat <- as.matrix(bdiag(list(matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25))))
#' res <- preprocess.genetic.data(case[, 1:10], father.genetic.data = dad[ , 1:10],
#'                                mother.genetic.data = mom[ , 1:10],
#'                                block.ld.mat = block.ld.mat[ , 1:10])
#'
#' @importFrom matrixStats colSds rowMaxs
#' @importFrom data.table data.table rbindlist setorder
#' @importFrom stats rbinom sd
#' @importFrom survival clogit
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom survival clogit strata coxph Surv
#' @export

preprocess.genetic.data <- function(case.genetic.data, complement.genetic.data = NULL, father.genetic.data = NULL,
    mother.genetic.data = NULL, block.ld.mat, min.allele.freq = 0, bp.param = bpparam(), snp.sampling.probs = NULL) {

    # make sure the appropriate genetic data is included
    if (is.null(complement.genetic.data) & (is.null(father.genetic.data) | is.null(mother.genetic.data))) {

        stop("Must include complement.genetic.data or both father.genetic.data and mother.genetic.data")

    }

    # make sure genotypes are integers
    if (!all(round(case.genetic.data) == case.genetic.data, na.rm = TRUE)){

        stop("case.genetic.data genotypes must be integers, not dosages imputed with uncertainty")

    }

    # check for missing data in the cases
    missing.case.geno <- is.na(case.genetic.data)

    ### find the snps with MAF < minimum threshold in the cases ###
    if (!is.null(father.genetic.data) & !is.null(mother.genetic.data)) {

        if (!all(round(father.genetic.data) == father.genetic.data, na.rm = TRUE)){

            stop("father.genetic.data genotypes must be integers, not dosages imputed with uncertainty")

        }

        if (!all(round(mother.genetic.data) == mother.genetic.data, na.rm = TRUE)){

            stop("mother.genetic.data genotypes must be integers, not dosages imputed with uncertainty")

        }

        # check for missing genotypes in moms and dads
        missing.father.geno <- is.na(father.genetic.data)
        missing.mother.geno <- is.na(mother.genetic.data)
        any.missing.geno <- missing.case.geno | missing.father.geno | missing.mother.geno

        # for now, set missing values to NA for the whole family
        case.genetic.data[any.missing.geno] <- NA
        father.genetic.data[any.missing.geno] <- NA
        mother.genetic.data[any.missing.geno] <- NA

        alt.allele.freqs <- colSums(father.genetic.data + mother.genetic.data, na.rm = TRUE)/(4 * colSums(!any.missing.geno))
        minor.alleles <- alt.allele.freqs < 0.5
        below.maf.threshold <- alt.allele.freqs > (1 - min.allele.freq) | alt.allele.freqs < min.allele.freq
        original.col.numbers <- which(!below.maf.threshold)
        names(original.col.numbers) <- NULL

        ### Compute the complement data ###
        complement.genetic.data <- father.genetic.data + mother.genetic.data - case.genetic.data

        ### recode the case and complement data so that 1 indicates a copy of the minor allele ###
        case.genetic.data[, !minor.alleles] <- 2 - case.genetic.data[, !minor.alleles]
        complement.genetic.data[, !minor.alleles] <- 2 - complement.genetic.data[, !minor.alleles]

        ### remove the snps not meeting the required allele frequency threshold ###
        case.genetic.data <- case.genetic.data[, !below.maf.threshold]
        complement.genetic.data <- complement.genetic.data[, !below.maf.threshold]
        block.ld.mat <- block.ld.mat[!below.maf.threshold, !below.maf.threshold]
        any.missing.geno <- any.missing.geno[ , !below.maf.threshold]

    } else if (!is.null(complement.genetic.data)) {

        if (!all(round(complement.genetic.data) == complement.genetic.data, na.rm = TRUE)){

            stop("complement.genetic.data genotypes must be integers, not dosages imputed with uncertainty")

        }

        # check for missing genotypes in complements
        missing.complement.geno <- is.na(complement.genetic.data)

        # for now, set any missing genotypes to NA for the family
        any.missing.geno <- missing.case.geno | missing.complement.geno
        case.genetic.data[any.missing.geno] <- NA
        complement.genetic.data[any.missing.geno] <- NA

        # identify minor alleles
        alt.allele.freqs <- colSums(case.genetic.data + complement.genetic.data, na.rm = TRUE)/(4 * colSums(!any.missing.geno))
        minor.alleles <- alt.allele.freqs < 0.5
        below.maf.threshold <- alt.allele.freqs > (1 - min.allele.freq) | alt.allele.freqs < min.allele.freq
        original.col.numbers <- which(!below.maf.threshold)
        names(original.col.numbers) <- NULL

        ### recode the case and complement data so that 1 indicates a copy of the minor allele ###
        case.genetic.data[, !minor.alleles] <- 2 - case.genetic.data[, !minor.alleles]
        complement.genetic.data[, !minor.alleles] <- 2 - complement.genetic.data[, !minor.alleles]

        ### remove the snps not meeting the required allele frequency threshold ###
        case.genetic.data <- case.genetic.data[, !below.maf.threshold]
        complement.genetic.data <- complement.genetic.data[, !below.maf.threshold]
        block.ld.mat <- block.ld.mat[!below.maf.threshold, !below.maf.threshold]
        any.missing.geno <- any.missing.geno[ , !below.maf.threshold]

    }

    if (is.null(snp.sampling.probs)){

        ### use conditional logistic regression to estimate univariate association ###
        case.status <- c(rep(1, nrow(case.genetic.data)), rep(0, nrow(complement.genetic.data)))
        ids <- rep(seq_len(nrow(case.genetic.data)), 2)

        res.list <- bplapply(seq_len(ncol(case.genetic.data)), function(snp, case.genetic.data, complement.genetic.data) {

            case.snp <- case.genetic.data[, snp]
            comp.snp <- complement.genetic.data[, snp]

            # get p-value of association from conditional logistic regression
            case.comp.geno <- c(case.snp, comp.snp)
            clogit.res <- clogit(case.status ~ case.comp.geno + strata(ids), method = "approximate")
            clogit.chisq <- summary(clogit.res)$logtest[1]

            return(list(case.snp = case.snp, comp.snp = comp.snp, chisq = clogit.chisq))

        }, case.genetic.data = case.genetic.data, complement.genetic.data = complement.genetic.data, BPPARAM = bp.param)
        chisq.stats <- do.call("c", lapply(res.list, function(x) x$chisq))

    } else {

        chisq.stats <- snp.sampling.probs[!below.maf.threshold]

    }

    # make sure case and comp data are stored correctly
    case.genetic.data <- as.matrix(case.genetic.data)
    storage.mode(case.genetic.data) <- "integer"
    complement.genetic.data <- as.matrix(complement.genetic.data)
    storage.mode(complement.genetic.data) <- "integer"

    # set any missing values to -9
    case.genetic.data[any.missing.geno] <- -9
    complement.genetic.data[any.missing.geno] <- -9

    return(list(case.genetic.data = case.genetic.data, complement.genetic.data = complement.genetic.data,
        chisq.stats = chisq.stats, original.col.numbers = original.col.numbers, block.ld.mat = block.ld.mat,
        minor.allele.vec = minor.alleles))

}
