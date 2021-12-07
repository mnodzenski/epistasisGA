#' A function to pre-process case-parent triad of affected/unaffected sibling data.
#'
#' This function performs several pre-processing steps, intended for use before function run.gadgets.
#'
#' @param case.genetic.data The genetic data of the disease affected children from case-parent trios or affected/unaffected sibling pairs.
#' Columns are SNP allele counts, and rows are individuals. This object may either be of class 'matrix' OR of class 'big.matrix'. If of class
#' 'big.matrix' it must be file backed as type 'integer' (see the bigmemory package for more information). The ordering of the columns must be consistent
#' with the LD structure specified in \code{ld.block.vec}. The genotypes cannot be dosages imputed with uncertainty.
#' @param complement.genetic.data A genetic dataset from the complements of the cases, where
#' \code{complement.genetic.data} = mother SNP counts + father SNP counts - case SNP counts. If using affected/unaffected siblings
#' this argument should be the genotypes for the unaffected siblings. This object may either be of class 'matrix' OR of class 'big.matrix'. If of class
#' 'big.matrix' it must be file backed as type 'integer' (see the bigmemory package for more information). Columns are SNP allele counts, rows are
#' families. If not specified, \code{father.genetic.data} and \code{mother.genetic.data} must be specified.
#' The genotypes cannot be dosages imputed with uncertainty.
#' @param father.genetic.data The genetic data for the fathers of the cases. Columns are SNP allele counts, rows are individuals.
#' Does not need to be specified if \code{complement.genetic.data} is specified. This object may either be of class 'matrix' OR of class 'big.matrix'. If of class
#' 'big.matrix' it must be file backed as type 'integer' (see the bigmemory package for more information). The genotypes cannot be dosages imputed with
#' uncertainty.
#' @param mother.genetic.data The genetic data for the mothers of the cases. Columns are SNP allele counts, rows are individuals.
#' Does not need to be specified if \code{complement.genetic.data} is specified. This object may either be of class 'matrix' OR of class 'big.matrix'. If of class
#' 'big.matrix' it must be file backed as type 'integer' (see the bigmemory package for more information). The genotypes cannot be dosages imputed with
#' uncertainty.
#' @param ld.block.vec An integer vector specifying the linkage blocks of the input SNPs. As an example, for 100 candidate SNPs, suppose
#' we specify \code{ld.block.vec <- c(25, 50, 25)}. This vector indicates that the input genetic data has 3 distinct linkage blocks, with
#' SNPs 1-25 in the first linkage block, 26-75 in the second block, and 76-100 in the third block. Note that this means the ordering of the columns (SNPs)
#' in \code{case.genetic.data} must be consistent with the LD blocks specified in \code{ld.block.vec}. In the absence of outside information,
#' a reasonable default is to consider SNPs to be in LD if they are located on the same biological chromosome. If not specified, this defaults
#' to assuming all input SNPs are in linkage, which may be overly conservative and could adversely affect performance.
#' @param bp.param The BPPARAM argument to be passed to bplapply when estimating marginal disease associations for each SNP.
#'  If using a cluster computer, this parameter needs to be set with care. See \code{BiocParallel::bplapply} for more details
#' @param snp.sampling.probs A vector indicating the sampling probabilities of the SNPs in \code{case.genetic.data}. SNPs will be sampled in the
#' genetic algorithm proportional to the values specified. If not specified, by default, chi-square statistics of association will be computed for
#' each SNP, and sampling will be proportional to the root of these statistics. If user specified, the  vector values need not sum to 1, they just need to be positive
#' real numbers. See argument \code{prob} from function \code{sample} for more details.
#' @param categorical.exposures A vector of integers corresponding to categorical exposures for the cases. Defaults to NULL,
#' which will result in GADGETS looking for epistatic interactions, rather than SNP by exposure interactions. \code{categorical.exposures}
#' should not be missing any data, families with missing exposure data should be removed from the analysis prior to input.
#' @param use.parents A integer indicating whether family level informativeness should be used alongside transmissions in computing GxE fitness scores. Defaults to 1,
#' indicating family level informativeness will be used. Specify 0 to only use transmission data.
#' @return A list containing the following:
#' \describe{
#'  \item{genetic.data.list}{A list of big.matrix.descriptor objects describing the locations of the input big.matrix objects
#'  containing the genetic data to be analyzed.}
#'  \item{chisq.stats}{A vector of chi-square statistics corresponding to marginal SNP-disease associations, if \code{snp.sampling.probs}
#'  is not specified, and \code{snp.sampling.probs} if specified.}
#'  \item{ld.block.vec}{A vector eaul to \code{cumsum(ld.block.vec)}.}
#'  \item{exposure}{A vector of categorical exposures, if specified, otherwise NULL.}
#' }
#'
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' case <- as.matrix(case)
#' dad <- as.matrix(dad)
#' mom <- as.matrix(mom)
#' res <- preprocess.genetic.data(case[, 1:10], father.genetic.data = dad[ , 1:10],
#'                                mother.genetic.data = mom[ , 1:10],
#'                                ld.block.vec = c(10))
#'
#' @importFrom bigmemory as.big.matrix describe attach.big.matrix
#' @importFrom data.table data.table rbindlist setorder
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom survival clogit strata coxph Surv
#' @export

preprocess.genetic.data <- function(case.genetic.data, complement.genetic.data = NULL, father.genetic.data = NULL,
    mother.genetic.data = NULL, ld.block.vec = NULL, bp.param = bpparam(), snp.sampling.probs = NULL,
    categorical.exposures = NULL, use.parents = 1) {

    #make sure the ld.block.vec is correctly specified
    if (is.null(ld.block.vec)){

        ld.block.vec <- ncol(case.genetic.data)

    } else {

        if (sum(ld.block.vec) != ncol(case.genetic.data)){

            stop("sum(ld.block.vec) must be equal to ncol(case.genetic.data)")

        }

    }

    # make sure the appropriate genetic data is included
    if (is.null(complement.genetic.data) & (is.null(father.genetic.data) | is.null(mother.genetic.data))) {

        stop("Must include complement.genetic.data or both father.genetic.data and mother.genetic.data")

    }

    ### if environmental exposures are provided, check their formatting ###
    if (!is.null(categorical.exposures)){

        # make sure categorical exposures are integers
        if (class(categorical.exposures) != "integer"){

            stop("categorical.exposures must be of class integer")

        }
        # identify families with missing exposure data
        missing.exposure <- is.na(categorical.exposures)
        if (sum(missing.exposure) > 0){

            stop(paste("Please remove", sum(missing.exposure), "families from analysis due to missing exposure(s)"))

        }

        # shorten the name of the exposures variable
        exposure <- categorical.exposures
        storage.mode(exposure) <- "integer"

        # get rid of any levels with only one case
        cases.per.level <- vapply(unique(exposure), function(exp.level){

            these.rows <- exposure == exp.level
            sum(these.rows)

        }, 1)

        one.case.levels <- unique(exposure)[cases.per.level == 1]
        if (length(one.case.levels) > 0){

            stop("At least two cases are required for each level of categorical.exposures.")

        }

        if (length(exposure) == 1){

            stop("exposure must have at least two levels")

        }

    } else {

        exposure <- NULL

    }

    # check formatting of input data and, if necessary, create memory mapped files
    if (!any(class(case.genetic.data) %in% c("matrix", "big.matrix"))){

        stop("case.genetic.data must be of class matrix or big.matrix")

    }

    if (any(class(case.genetic.data) == "matrix")){

        if (!all(round(case.genetic.data) == case.genetic.data, na.rm = TRUE)){

            stop("case.genetic.data genotypes must be integers, not dosages imputed with uncertainty")

        }

        storage.mode(case.genetic.data) <- "integer"

        # convert to big.matrix
        dimnames(case.genetic.data) <- NULL
        case.bm <- as.big.matrix(case.genetic.data, type = "integer")
        rm(case.genetic.data)

    } else if (class(case.genetic.data) == "big.matrix"){

        if (! describe(case.genetic.data)@description$type %in% c("integer")){

            stop("case.genetic.data must be a big.matrix of type integer. To convert, see function deepcopy from package bigmemory.")

        }

        if (describe(case.genetic.data)@description$sharedType != "FileBacked"){

            stop("case.genetic.data must be a file backed big.matrix (case.genetic.data@description$sharedType == 'FileBacked')")

        }

        case.bm <- case.genetic.data

    }

    if (!is.null(complement.genetic.data) & !any(class(complement.genetic.data) %in% c("matrix", "big.matrix"))){

        stop("complement.genetic.data must be of class matrix or big.matrix")

    }

    if (!is.null(complement.genetic.data) & any(class(complement.genetic.data) == "matrix")){

        if (!all(round(complement.genetic.data) == complement.genetic.data, na.rm = TRUE)){

            stop("complement.genetic.data genotypes must be integers, not dosages imputed with uncertainty")

        }

        storage.mode(complement.genetic.data) <- "integer"

        # convert to big.matrix
        dimnames(complement.genetic.data) <- NULL
        comp.bm <- as.big.matrix(complement.genetic.data, type = "integer")
        rm(complement.genetic.data)

    } else if (!is.null(complement.genetic.data) & class(complement.genetic.data) == "big.matrix"){

        if (! describe(complement.genetic.data)@description$type %in% c("integer")){

            stop("complement.genetic.data must be a big.matrix of type integer. To convert, see function deepcopy from package bigmemory.")

        }

        if (describe(complement.genetic.data)@description$sharedType != "FileBacked"){

            stop("complement.genetic.data must be a file backed big.matrix (complement.genetic.data@description$sharedType == 'FileBacked')")

        }

        comp.bm <- complement.genetic.data

    }

    if (!is.null(mother.genetic.data) & !any(class(mother.genetic.data) %in% c("matrix", "big.matrix"))){

        stop("mother.genetic.data must be of class matrix or big.matrix")

    }

    if (!is.null(mother.genetic.data) & any(class(mother.genetic.data) == "matrix")){

        if (!all(round(mother.genetic.data) == mother.genetic.data, na.rm = TRUE)){

            stop("mother.genetic.data genotypes must be integers, not dosages imputed with uncertainty")

        }

        storage.mode(mother.genetic.data) <- "integer"

        # convert to big.matrix
        dimnames(mother.genetic.data) <- NULL
        mother.bm <- as.big.matrix(mother.genetic.data, type = "integer")
        rm(mother.genetic.data)

    } else if (!is.null(mother.genetic.data) & class(mother.genetic.data) == "big.matrix"){

        if (! describe(mother.genetic.data)@description$type %in% c("integer")){

            stop("mother.genetic.data must be a big.matrix of type integer. To convert, see function deepcopy from package bigmemory.")

        }

        if (describe(mother.genetic.data)@description$sharedType != "FileBacked"){

            stop("mother.genetic.data must be a file backed big.matrix (mother.genetic.data@description$sharedType == 'FileBacked')")

        }

        mother.bm <- mother.genetic.data

    }

    if (!is.null(father.genetic.data) & !any(class(father.genetic.data) %in% c("matrix", "big.matrix"))){

        stop("father.genetic.data must be of class matrix or big.matrix")

    }

    if (!is.null(father.genetic.data) & any(class(father.genetic.data) == "matrix")){

        if (!all(round(father.genetic.data) == father.genetic.data, na.rm = TRUE)){

            stop("father.genetic.data genotypes must be integers, not dosages imputed with uncertainty")

        }
        storage.mode(father.genetic.data) <- "integer"

        # convert to big.matrix
        dimnames(father.genetic.data) <- NULL
        father.bm <- as.big.matrix(father.genetic.data, type = "integer")
        rm(father.genetic.data)

    } else if (!is.null(father.genetic.data) & class(father.genetic.data) == "big.matrix"){

        if (! describe(father.genetic.data)@description$type %in% c("integer")){

            stop("father.genetic.data must be a big.matrix of type integer. To convert, see function deepcopy from package bigmemory.")

        }

        if (describe(father.genetic.data)@description$sharedType != "FileBacked"){

            stop("father.genetic.data must be a file backed big.matrix (father.genetic.data@description$sharedType == 'FileBacked')")

        }

        father.bm <- father.genetic.data

    }

    # make a list of big matrix objects
    if (exists("comp.bm")){

        bm.list <- list(case = case.bm, complement = comp.bm)

    } else {

        bm.list <- list(case = case.bm, mother = mother.bm, father = father.bm)

    }

    # if needed compute sampling probs
    if (is.null(snp.sampling.probs)){

        ### use conditional logistic regression to estimate univariate association ###
        n.fam <- nrow(case.bm)
        n.candidate.snps <- ncol(case.bm)
        case.status <- c(rep(1, n.fam), rep(0, n.fam))
        ids <- rep(seq_len(n.fam), 2)

        if (is.null(categorical.exposures)){

            res.list <- bplapply(seq_len(n.candidate.snps), function(snp, bm.list) {

                case.snp <- bm.list$case[ , snp]
                if (length(bm.list) == 3){

                    mom.snp <- bm.list$mother[ , snp]
                    dad.snp <- bm.list$father[ , snp]
                    comp.snp <- mom.snp + dad.snp - case.snp

                } else {

                    comp.snp <- bm.list$complement[ , snp]

                }

                # get p-value of association from conditional logistic regression
                case.comp.geno <- c(case.snp, comp.snp)
                clogit.res <- clogit(case.status ~ case.comp.geno + strata(ids), method = "approximate")
                clogit.chisq <- summary(clogit.res)$logtest[1]

                return(list(case.snp = case.snp, comp.snp = comp.snp, chisq = clogit.chisq))

            }, bm.list = bm.list, BPPARAM = bp.param)
            chisq.stats <- do.call("c", lapply(res.list, function(x) x$chisq))

        } else {

            if (use.parents != 2){

                exposure.var <- factor(rep(exposure, 2))
                res.list <- bplapply(seq_len(n.candidate.snps), function(snp, bm.list, exposure.var) {

                    case.snp <- bm.list$case[ , snp]
                    if (length(bm.list) == 3){

                        mom.snp <- bm.list$mother[ , snp]
                        dad.snp <- bm.list$father[ , snp]
                        comp.snp <- mom.snp + dad.snp - case.snp

                    } else {

                        comp.snp <- bm.list$complement[ , snp]

                    }

                    # get p-value of snp-exposure association from conditional logistic regression
                    case.comp.geno <- c(case.snp, comp.snp)
                    df <- data.table(case.status = case.status, case.comp.geno = case.comp.geno, exposure = exposure.var, ids = ids)
                    full.model <- clogit(case.status ~ case.comp.geno + case.comp.geno:exposure + strata(ids), method = "approximate", data  = df)
                    full.model.ll <- full.model$loglik[2]
                    reduced.model <- clogit(case.status ~ case.comp.geno + strata(ids), method = "approximate", data  = df)
                    reduced.model.ll <- reduced.model$loglik[2]
                    clogit.chisq <- 2*(full.model.ll - reduced.model.ll)

                    return(list(case.snp = case.snp, comp.snp = comp.snp, chisq = clogit.chisq))

                }, bm.list = bm.list, exposure.var = exposure.var, BPPARAM = bp.param)
                chisq.stats <- do.call("c", lapply(res.list, function(x) x$chisq))

            } else {

                exposure.var <- as.factor(exposure)
                res.list <- bplapply(seq_len(n.candidate.snps), function(snp, bm.list, exposure.var) {

                    case.snp <- bm.list$case[ , snp]
                    if (length(bm.list) == 3){

                        mom.snp <- bm.list$mother[ , snp]
                        dad.snp <- bm.list$father[ , snp]
                        comp.snp <- mom.snp + dad.snp - case.snp

                    } else {

                        comp.snp <- bm.list$complement[ , snp]

                    }

                    # get p-value of snp-exposure association from conditional logistic regression
                    case.comp.diff <- as.numeric(case.snp != comp.snp)
                    glm.res <- glm(case.comp.diff ~ exposure.var, family = "binomial")
                    test.stat <- summary(glm.res)$null.deviance - summary(glm.res)$deviance

                    return(list(chisq = test.stat))

                }, bm.list = bm.list, exposure.var = exposure.var, BPPARAM = bp.param)
                chisq.stats <- do.call("c", lapply(res.list, function(x) x$chisq))

            }

        }

    }

    # take cumulative sum of ld.block.vec for output
    out.ld.vec <- cumsum(ld.block.vec)
    storage.mode(out.ld.vec) <- "integer"

    #### clean up chisq stats for models that did not converge ###
    chisq.stats[chisq.stats <= 0] <- 10^-10
    chisq.stats[is.infinite(chisq.stats)] <- max(chisq.stats[is.finite(chisq.stats)])

    ### if running GxE create required inputs ###
    if (!is.null(exposure)){

        exposure.levels <- unique(exposure)
        storage.mode(exposure.levels) <- "integer"

    } else {

        exposure.levels <- NULL

    }

    if (!"complement" %in% names(bm.list)){

        comp.data <- mother.bm[] + father.bm[] - case.bm[]

    } else {

        comp.data <- comp.bm[]

    }

    case.data <- case.bm[]

    # set missing to -9
    if (any(is.na(case.data)) | any(is.na(comp.data))){

        case.data[is.na(case.data) | is.na(comp.data)] <- -9
        comp.data[is.na(case.data) | is.na(comp.data)] <- -9

    }

    return(list(case.genetic.data = case.data, complement.genetic.data = comp.data, chisq.stats = chisq.stats, ld.block.vec = out.ld.vec,
        exposure = exposure, exposure.levels = exposure.levels, use.parents = use.parents))

}
