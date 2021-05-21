#' A function to pre-process case-parent triad of affected/unaffected sibling data.
#'
#' This function performs several pre-processing steps, intended for use before function run.gadgets.
#'
#' @param case.genetic.data.trio.study The genetic data of the disease affected children from case-parent trios. Columns are SNP allele counts, and rows are individuals.
#' The ordering of the columns must be consistent with the LD structure specified in \code{block.ld.mat}. If using a mix of trio and sibling data, specify the trio
#' data here and in \code{father.genetic.data} and \code{mother.genetic.data}. Specify the the sibling study data in argument \code{case.genetic.data.sibling.study} and
#' \code{sibling.genetic.data}.
#' @param father.genetic.data The genetic data for the fathers of the cases in \code{case.genetic.data.trio.study}. Columns are SNP allele counts, rows are individuals.
#' @param mother.genetic.data The genetic data for the mothers of the cases in \code{case.genetic.data.trio.study}. Columns are SNP allele counts, rows are individuals.
#' @param complement.genetic.data.trio.study A genetic dataset from the complements of the cases in \code{case.genetic.data.trio.study}, where
#' \code{complement.genetic.data.trio.study} = mother SNP counts + father SNP counts - case SNP counts.
#' Columns are SNP allele counts, rows are families. Alternatively, \code{father.genetic.data} and \code{mother.genetic.data} can be specified.
#' @param case.covars.trio.study A matrix or data.frame of binary covariates to be included in the stochastic search corresponding to the affected cases in
#' \code{case.genetic.data.trio.study}. Columns correspond to different covariates, rows correspond to different families.
#' @param case.genetic.data.sibling.study The genetic data of the disease affected children from affected/unaffected sibling pairs. Columns are SNP allele counts, and rows are individuals.
#' The ordering of the columns must be consistent with the LD structure specified in \code{block.ld.mat}.
#' @param sibling.genetic.data The genetic data of the disease affected children from affected/unaffected sibling pairs. Columns are SNP allele counts, and rows are individuals.
#' The ordering of the columns must be consistent with the LD structure specified in \code{block.ld.mat}. If desired, for case-parent trio studies,
#' this argument can be entered as \code{sibling.genetic.data} = mother SNP counts + father SNP counts - case SNP counts. This would be needed, for example, if
#' any of the input SNPs are on the X-chromosome. Columns are SNP allele counts, rows are families.
#' @param case.covars.sibling.study A matrix or data.frame of binary covariates to be included in the stochastic search corresponding to the affected cases in
#' \code{case.genetic.data.sibling.study}. Columns correspond to different covariates, rows correspond to different families.
#' @param sibling.covars A matrix or data.frame of sibling covariates, with columns matching those specified in \code{case.covars.sibling.study}.
#' @param block.ld.mat A logical, block diagonal matrix indicating whether the SNPs in \code{case.genetic.data} should be considered
#'  to be in linkage disequilibrium. Note that this means the ordering of the columns (SNPs) in \code{case.genetic.data} must be consistent
#'  with the LD blocks specified in \code{ld.block.mat}. In the absence of outside information, a reasonable default is to consider SNPs
#'  to be in LD if they are located on the same biological chromosome. If not specified, this defaults to the assumption that all SNPs are potentially in linkage,
#'  which may adversely affect performance compared to a more carefully and biologically relevant specification.
#' @param min.allele.freq The minimum minor allele frequency required for a SNP to be considered for inclusion in the genetic algorithm.
#' Any SNPs with MAF < \code{min.allele.freq} in the parents, or the combined group of affected and unaffected siblings, will be filtered out. Defaults to 0 (no filtering).
#' @param bp.param The BPPARAM argument to be passed to bplapply when estimating marginal disease associations for each SNP.
#'  If using a cluster computer, this parameter needs to be set with care. See \code{BiocParallel::bplapply} for more details
#' @param snp.sampling.probs A vector indicating the sampling probabilities of the SNPs in \code{case.genetic.data}. SNPs will be sampled in the
#' genetic algorithm proportional to the values specified. If not specified, by default, chi-square statistics of association will be computed for
#' each SNP, and sampling will be proportional to the root of these statistics. If covariates are included, the SNPs will be sampled proportional to the
#' chi-squared statistic corresponding to the SNP-covariate interaction. If user specified, the  vector values need not sum to 1, they just need to be positive
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
#'  \item{exposure}{A vector of categorical exposures, if specified, otherwise NULL.}
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

preprocess.genetic.data <- function(case.genetic.data.trio.study = NULL, father.genetic.data = NULL,
    mother.genetic.data = NULL, complement.genetic.data.trio.study = NULL, case.covars.trio.study = NULL,
    case.genetic.data.sibling.study = NULL, sibling.genetic.data = NULL, case.covars.sibling.study = NULL,
    sibling.covars = NULL, block.ld.mat = NULL, min.allele.freq = 0, bp.param = bpparam(), snp.sampling.probs = NULL) {


    # make sure the appropriate genetic data is included
    if (is.null(case.genetic.data.trio.study) & is.null(case.genetic.data.sibling.study)) {

        stop("Must include either case.genetic.data.trio.study or case.genetic.data.sibling.study")

    }

    if (!is.null(case.genetic.data.trio.study) & ((is.null(father.genetic.data) | is.null(mother.genetic.data)) |
                                                  is.null(complement.genetic.data.trio.study))) {

        stop("Must include father.genetic.data and mother.genetic.data or complement.genetic.data.trio.study corresponding to case.genetic.data.trio.study")

    }

    if (!is.null(case.genetic.data.sibling.study) & (is.null(sibling.genetic.data))) {

        stop("Must include sibling.genetic.data corresponding to case.genetic.data.sibling.study")

    }

    if (!is.null(case.genetic.data.trio.study) & !is.null(case.genetic.data.sibling.study)){

        if (ncol(case.genetic.data.sibling.study) != ncol(case.genetic.data.trio.study)){

            stop("case.genetic.data.sibling.study and case.genetic.data.trio.study do not contain the same number of SNPs")

        }

    }

    # check covariates
    if (!is.null(sibling.covars) & (is.null(case.covars.sibling.study))) {

        stop("Must include case.covars.sibling.study corresponding to sibling.covars")

    }

    if (is.null(sibling.covars) & (!is.null(case.covars.sibling.study))) {

        stop("Must include sibling.covars corresponding to case.covars.sibling.study")

    }

    if (!is.null(case.covars.trio.study) & !is.null(case.covars.sibling.study)){

        if (ncol(case.covars.trio.study) != ncol(case.covars.sibling.study)){

            stop("case.covars.trio.study and case.covars.sibling.study do not contain the same number of covariates")

        }

    }

    ### if covariates are provided, make sure they are binary ###
    use.covars <- FALSE
    if (!is.null(case.covars.trio.study)){

        use.covars <- TRUE
        if (any(!case.covars.trio.study %in% c(0, 1))){

            stop("case.covars.trio.study must be coded 0 or 1")

        }

        # check for missing case covariates
        missing.case.covars.trio <- is.na(case.covars.trio.study)

        # construct the complement covars
        comp.covars <- 1 - case.covars.trio.study

        # temporarily set covariate missing for either case or complement where appropriate
        case.genetic.data.trio.study[missing.case.covars.trio] <- NA
        comp.covars[missing.case.covars.trio] <- NA

    }

    if (!is.null(case.covars.sibling.study)){

        use.covars <- TRUE
        if (any(!case.covars.sibling.study %in% c(0, 1, NA))){

            stop("case.covars.sibling.study must be coded 0 or 1 or NA")

        }

        if (any(!sibling.covars %in% c(0, 1, NA))){

            stop("sibling.covars must be coded 0 or 1 or NA")

        }

        # check for missing case covariates
        missing.case.covars.sibling.study <- is.na(case.covars.sibling.study)

        # missing sibling covariates
        missing.sibling.covars <- is.na(sibling.covars)
        either.missing.covars <- missing.case.covars.sibling.study | missing.sibling.covars

        # set missing where either case or sibling is missing covars (i.e., uninformative)
        case.covars.sibling.study[either.missing.covars] <- NA
        sibling.covars[either.missing.covars] <- NA

    }

    # initiate block.ld.mat if not specified
    if (is.null(block.ld.mat)){

        if (!is.null(case.genetic.data.trio.study)){

            block.ld.mat <- matrix(TRUE, ncol(case.genetic.data.trio.study), ncol(case.genetic.data.trio.study))

        } else if (!is.null(case.genetic.data.sibling.study)) {

            block.ld.mat <- matrix(TRUE, ncol(case.genetic.data.sibling.study), ncol(case.genetic.data.sibling.study))

        }


    }

    # combine into case and complement covar data frames
    # note we need a different version of the covars for the LRT, where the complement
    # gets the same covars as the case (by convention of those models)
    if (!is.null(case.covars.trio.study) & !is.null(case.covars.sibling.study)){

        case.covars <- rbind(case.covars.trio.study, case.covars.sibling.study)
        comp.covars <- rbind(comp.covars, sibling.covars)
        comp.covars.model <- rbind(case.covars.trio.study, sibling.covars)
        any.missing.covars <- rbind(missing.case.covars.trio, either.missing.covars)

    } else if (!is.null(case.covars.trio.study) & is.null(case.covars.sibling.study)){

        case.covars <- case.covars.trio.study
        comp.covars.model <- case.covars.trio.study
        any.missing.covars <- missing.case.covars.trio

    } else if (is.null(case.covars.trio.study) & !is.null(case.covars.sibling.study)){

        case.covars <- case.covars.sibling.study
        comp.covars <- sibling.covars
        comp.covars.model <- sibling.covars
        any.missing.covars <- either.missing.covars

    }

    ### process trio genetic data ###
    if (!is.null(case.genetic.data.trio.study)) {

        # make sure genotypes are integers
        if (!all(round(case.genetic.data.trio.study) == case.genetic.data.trio.study, na.rm = TRUE)){

            stop("case.genetic.data.trio.study genotypes must be integers, not dosages imputed with uncertainty")

        }

        # check for missing genotypes
        missing.case.trio.geno <- is.na(case.genetic.data.trio.study)

        if (!is.null(father.genetic.data) & !is.null(mother.genetic.data)){

            if (!all(round(father.genetic.data) == father.genetic.data, na.rm = TRUE)){

                stop("father.genetic.data genotypes must be integers, not dosages imputed with uncertainty")

            }

            if (!all(round(mother.genetic.data) == mother.genetic.data, na.rm = TRUE)){

                stop("mother.genetic.data genotypes must be integers, not dosages imputed with uncertainty")

            }

            missing.father.geno <- is.na(father.genetic.data)
            missing.mother.geno <- is.na(mother.genetic.data)
            any.missing.geno <- missing.case.trio.geno | missing.father.geno | missing.mother.geno

            # for now, set missing values to NA for the whole family
            case.genetic.data.trio.study[any.missing.geno] <- NA
            father.genetic.data[any.missing.geno] <- NA
            mother.genetic.data[any.missing.geno] <- NA

            ### Compute the complement data ###
            complement.genetic.data <- father.genetic.data + mother.genetic.data - case.genetic.data

        } else if (!is.null(complement.genetic.data.trio.study)){

            complement.genetic.data <- complement.genetic.data.trio.study

            if (!all(round(complement.genetic.data.trio.study) == complement.genetic.data.trio.study, na.rm = TRUE)){

                stop("complement.genetic.data.trio.study genotypes must be integers, not dosages imputed with uncertainty")

            }

            missing.complement.geno <- is.na(complement.genetic.data.trio.study)
            any.missing.geno <- missing.case.trio.geno | missing.complement.geno

            # for now, set missing values to NA for the whole family
            case.genetic.data.trio.study[any.missing.geno] <- NA
            complement.genetic.data[any.missing.geno] <- NA
        }

        if (is.null(case.genetic.data.sibling.study)){

            case.genetic.data <- case.genetic.data.trio.study

            ### only for trio studies, if mixed, need to combine data first), identify minor alleles ###
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

            ### indicator that these are from trio study data ###
            trio.data <- rep(TRUE, nrow(case.genetic.data))

        }


    }

    if (!is.null(case.genetic.data.sibling.study)) {

        ### require integer genotypes ###
        if (!all(round(case.genetic.data.sibling.study) == case.genetic.data.sibling.study, na.rm = TRUE)){

            stop("case.genetic.data.sibling.study genotypes must be integers, not dosages imputed with uncertainty")

        }

        if (!all(round(sibling.genetic.data) == sibling.genetic.data, na.rm = TRUE)){

            stop("sibling.genetic.data genotypes must be integers, not dosages imputed with uncertainty")

        }

        # check for missing genotypes
        missing.case.geno.sibling.study <- is.na(case.genetic.data.sibling.study)
        missing.sibling.geno <- is.na(sibling.genetic.data)

        # for now, set any missing genotypes to NA for the family
        any.missing.geno <- missing.case.geno.sibling.study | missing.sibling.geno
        case.genetic.data.sibling.study[any.missing.geno] <- NA
        sibling.genetic.data[any.missing.geno] <- NA

        if (is.null(case.genetic.data.trio.study)){

            case.genetic.data <- case.genetic.data.sibling.study
            complement.genetic.data <- sibling.genetic.data

            ### indicator that these are sibling data ###
            trio.data <- rep(FALSE, nrow(case.genetic.data))

        } else {

            case.genetic.data <- rbind(case.genetic.data.trio.study, case.genetic.data.sibling.study)
            complement.genetic.data <- rbind(complement.genetic.data, sibling.genetic.data)
            any.missing.geno <- is.na(case.genetic.data) | is.na(complement.genetic.data)

            ### indicator that these are mixed trio/sibling data ###
            trio.data <- c(rep(TRUE, nrow(case.genetic.data.trio.study)),
                           rep(FALSE, nrow(case.genetic.data.sibling.study)))
        }

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

        if (!use.covars){

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

            res.list <- bplapply(seq_len(ncol(case.genetic.data)), function(snp, case.genetic.data, complement.genetic.data,
                                                                            case.covars, comp.covars.model) {

                case.snp <- case.genetic.data[, snp]
                comp.snp <- complement.genetic.data[, snp]

                # get p-value of snp-exposure association from conditional logistic regression
                case.comp.geno <- c(case.snp, comp.snp)
                df<- data.frame(rbind(case.covars, comp.covars.model))
                covar.cols <- paste0("covar", ncol(df))
                colnames(df) <- covar.cols
                df$case.status <- case.status
                df$case.comp.geno <- case.comp.geno
                df$ids <- ids

                # model specification

                # using only main effect of genotypes and any terms including gene/covar interaction
                inter.orders <- seq(2, length(covar.cols) + 1)
                inter.terms <- unlist(lapply(inter.orders, function(int.order){

                    order.int.terms <- combn(c("case.comp.geno", covar.cols), int.order)
                    int.terms <- apply(order.int.terms, 2, function(x) paste(x, collapse = "*"))
                    geno.int.terms <- grep("case.comp.geno", int.terms, value = TRUE)
                    return(geno.int.terms)

                }))
                inter.terms <- paste(inter.terms, collapse = "+")
                full.model.formula <- as.formula(paste("case.status ~ case.comp.geno", inter.terms, "strata(ids)", sep = "+"))
                full.model <- clogit(full.model.formula, method = "approximate", data  = df)
                full.model.ll <- full.model$loglik
                reduced.model <- clogit(case.status ~ case.comp.geno + strata(ids), method = "approximate", data  = df)
                reduced.model.ll <- reduced.model$loglik
                clogit.chisq <- 2*log(full.model.ll - reduced.model.ll)

                return(list(case.snp = case.snp, comp.snp = comp.snp, chisq = clogit.chisq))

            }, case.genetic.data = case.genetic.data, complement.genetic.data = complement.genetic.data, case.covars = case.covars,
            comp.covars.model = comp.covars.model, BPPARAM = bp.param)
            chisq.stats <- do.call("c", lapply(res.list, function(x) x$chisq))

        }

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

    if (use.covars){

        case.covars <- as.matrix(case.covars)
        storage.mode(case.covars) <- "integer"
        comp.covars <- as.matrix(comp.covars)
        storage.mode(comp.covars) <- "integer"
        case.covars[any.missing.covars] <- -9
        comp.covars[any.missing.covars] <- -9


    } else {

        case.covars <- matrix(1, 1, 1, drop = FALSE)
        storage.mode(case.covars) <- "integer"
        comp.covars <- matrix(1, 1, 1, drop = FALSE)

    }

    return(list(case.genetic.data = case.genetic.data, complement.genetic.data = complement.genetic.data,
                chisq.stats = chisq.stats, original.col.numbers = original.col.numbers, block.ld.mat = block.ld.mat,
                minor.allele.vec = minor.alleles, case.covars = case.covars, comp.covars = comp.covars,
                use.covars = use.covars, trio.data = trio.data))
}




