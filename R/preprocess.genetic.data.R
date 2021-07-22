#' A function to pre-process case-parent triad of affected/unaffected sibling data.
#'
#' This function performs several pre-processing steps, intended for use before function run.gadgets.
#'
#' @param case.genetic.data The genetic data of the disease affected children from case-parent trios or affected/unaffected sibling pairs.
#' Columns are SNP allele counts, and rows are individuals. This object may either be of class 'matrix' OR of class 'big.matrix'. If of class
#' 'big.matrix' it must be file backed as type 'short' (see the bigmemory package for more information). The ordering of the columns must be consistent
#' with the LD structure specified in \code{ld.block.vec}. The genotypes cannot be dosages imputed with uncertainty. If any data are missing for a particular
#' family for a particular SNP, that SNP's genotype should be coded as -9 for the entire family, (\code{case.genetic.data} and
#' \code{father.genetic.data}/\code{mother.genetic.data} or \code{case.genetic.data} and \code{complement.genetic.data}).
#' @param complement.genetic.data A genetic dataset from the complements of the cases, where
#' \code{complement.genetic.data} = mother SNP counts + father SNP counts - case SNP counts. If using affected/unaffected siblings
#' this argument should be the genotypes for the unaffected siblings. This object may either be of class 'matrix' OR of class 'big.matrix'. If of class
#' 'big.matrix' it must be file backed as type 'short' (see the bigmemory package for more information). Columns are SNP allele counts, rows are
#' families. If not specified, \code{father.genetic.data} and \code{mother.genetic.data} must be specified.
#' The genotypes cannot be dosages imputed with uncertainty. If any data are missing for a particular family for a particular SNP, that SNP's genotype
#' should be coded as -9 for the entire family (\code{case.genetic.data} and \code{complement.genetic.data}).
#' @param father.genetic.data The genetic data for the fathers of the cases. Columns are SNP allele counts, rows are individuals.
#' Does not need to be specified if \code{complement.genetic.data} is specified. This object may either be of class 'matrix' OR of class 'big.matrix'. If of class
#' 'big.matrix' it must be file backed as type 'short' (see the bigmemory package for more information). The genotypes cannot be dosages imputed with
#' uncertainty. If any data are missing for a particular family for a particular SNP, that SNP's genotype should be coded as -9 for the entire family,
#' (\code{case.genetic.data} and \code{father.genetic.data}/\code{mother.genetic.data}).
#' @param mother.genetic.data The genetic data for the mothers of the cases. Columns are SNP allele counts, rows are individuals.
#' Does not need to be specified if \code{complement.genetic.data} is specified. This object may either be of class 'matrix' OR of class 'big.matrix'. If of class
#' 'big.matrix' it must be file backed as type 'short' (see the bigmemory package for more information). The genotypes cannot be dosages imputed with
#' uncertainty. If any data are missing for a particular family for a particular SNP, that SNP's genotype should be coded as -9 for the entire family,
#' (\code{case.genetic.data} and \code{father.genetic.data}/\code{mother.genetic.data}).
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
#' @param categorical.exposures.risk.ranks An optional named list indicating the hypothesized relationship to risk
#' among the levels of \code{categorical.exposures}. The number of list elements must be equal to the number
#' of distinct levels of \code{categorical.exposures} and the list element names should be the
#' distinct levels of \code{categorical.exposures}. The list element values should be integers corresponding to
#' the relative rank of hypothesized risk corresponding to an exposure, with 1 corresponding to the lowest risk
#' level. For example, suppose \code{categorical.exposures} has two levels, 1 and 2, and an analyst is interested
#' only in identifying SNPs that are synergistically risk-related in the presence of exposure level 2. The analyst
#' should specify \code{list("1" = 1, "2" = 2)} for \code{categorical.exposures.risk.ranks}. Similarly, for
#' an exposure with levels 1, 2, and 3, with hypothesized increasing risk relevance with each level, an
#' analyst could specify \code{list("1" = 1, "2" = 2, "3" = 3)}. See the package vignette for more detailed
#' examples. If not specified, no risk-related ordering is assumed among the levels of \code{categorical.exposures}.
#' @param big.matrix.file.path  This argument specifies a directory where memory mapped files of class 'big.memory'
#' will be saved on disk for use in running the GADGETS method, allowing use of genetic datasets that do no
#' fit into RAM. This argument must be specified if (1) \code{case.genetic.data} or \code{complement.genetic.data}
#' is not a file backed big.matrix
#' (see package bigmemory) or (2) \code{mother.genetic.data} and \code{father.genetic.data} are specified
#' and \code{complement.genetic.data} is not specified.
#'
#' @return A list containing the following:
#' \describe{
#'  \item{genetic.data.list}{A list of big.matrix.descriptor objects describing the locations of the input big.matrix objects
#'  containing the genetic data to be analyzed.}
#'  \item{chisq.stats}{A vector of chi-square statistics corresponding to marginal SNP-disease associations, if \code{snp.sampling.probs}
#'  is not specified, and \code{snp.sampling.probs} if specified.}
#'  \item{ld.block.vec}{A vector eaul to \code{cumsum(ld.block.vec)}.}
#'  \item{exposure}{A vector of categorical exposures, if specified, otherwise NULL.}
#'  \item{exposure.risk.levels}{The list specified in input argument categorical.exposures.risk.ranks.}
#' }
#'
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' res <- preprocess.genetic.data(case[, 1:10], father.genetic.data = dad[ , 1:10],
#'                                mother.genetic.data = mom[ , 1:10],
#'                                ld.block.vec = c(10),
#'                                big.matrix.file.path = "tmp")
#' unlink(tmp)
#'
#' @importFrom bigmemory as.big.matrix describe attach.big.matrix deepcopy
#' @importFrom data.table data.table rbindlist setorder
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom survival clogit strata coxph Surv
#' @export

preprocess.genetic.data <- function(case.genetic.data, complement.genetic.data = NULL, father.genetic.data = NULL,
    mother.genetic.data = NULL, ld.block.vec = NULL, bp.param = bpparam(), snp.sampling.probs = NULL,
    categorical.exposures = NULL, categorical.exposures.risk.ranks = NULL, big.matrix.file.path = NULL) {

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

        # if specified, make sure the hypothesized risk vector contains the correct number
        # of elements
        if (!is.null(categorical.exposures.risk.ranks)){

            rr.names <- names(categorical.exposures.risk.ranks)
            unique.exposures <- unique(categorical.exposures)
            correct.rr <- all(unique.exposures %in% rr.names)
            if (!correct.rr){

                stop("names of categorical.exposures.risk.ranks must match the unique elements of categorical.exposures")

            }

        }

        # identify families with missing exposure data
        missing.exposure <- is.na(categorical.exposures)
        if (sum(missing.exposure) > 0){

            stop(paste("Please remove", sum(missing.exposure), "families from analysis due to missing exposure(s)"))

        }

        # shorten the name of the exposures variable
        exposure <- categorical.exposures
        exposure.risk.levels <- categorical.exposures.risk.ranks
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

        storage.mode(case.genetic.data) <- "integer"

        if (is.null(big.matrix.file.path)){

            stop("please specify big.matrix.file.path")

        }

        # convert to big.matrix
        if (!dir.exists(big.matrix.file.path)){

            dir.create(big.matrix.file.path, recursive = TRUE)

        }
        dimnames(case.genetic.data) <- NULL
        big.matrix.file.path <- normalizePath(big.matrix.file.path)
        case.bm <- as.big.matrix(case.genetic.data, type = "short", backingfile = "case_bm",
                                 backingpath = big.matrix.file.path, descriptorfile = "case_bm_desc.rds",
                                 binarydescriptor = TRUE)
        rm(case.genetic.data)

    } else if (class(case.genetic.data) == "big.matrix"){

        if (! describe(case.genetic.data)@description$type %in% c("short")){

            stop("case.genetic.data must be a big.matrix of type short. To convert, see function deepcopy from package bigmemory.")

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

        storage.mode(complement.genetic.data) <- "integer"

        if (is.null(big.matrix.file.path)){

            stop("please specify big.matrix.file.path")

        }

        # convert to big.matrix
        dimnames(complement.genetic.data) <- NULL
        comp.bm <- as.big.matrix(complement.genetic.data, type = "short", backingfile = "comp_bm",
                                 backingpath = big.matrix.file.path, descriptorfile = "comp_bm_desc.rds",
                                 binarydescriptor = TRUE)
        rm(complement.genetic.data)

    } else if (!is.null(complement.genetic.data) & class(complement.genetic.data) == "big.matrix"){

        if (! describe(complement.genetic.data)@description$type %in% c("short")){

            stop("complement.genetic.data must be a big.matrix of type short. To convert, see function deepcopy from package bigmemory.")

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

        storage.mode(mother.genetic.data) <- "integer"

        if (is.null(big.matrix.file.path)){

            stop("please specify big.matrix.file.path")

        }

        # convert to big.matrix
        dimnames(mother.genetic.data) <- NULL
        mother.bm <- as.big.matrix(mother.genetic.data, type = "short", backingfile = "mother_bm",
                                 backingpath = big.matrix.file.path, descriptorfile = "mother_bm_desc.rds",
                                 binarydescriptor = TRUE)
        rm(mother.genetic.data)

    } else if (!is.null(mother.genetic.data) & class(mother.genetic.data) == "big.matrix"){

        if (! describe(mother.genetic.data)@description$type %in% c("short")){

            stop("mother.genetic.data must be a big.matrix of type short. To convert, see function deepcopy from package bigmemory.")

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

        storage.mode(father.genetic.data) <- "integer"

        if (is.null(big.matrix.file.path)){

            stop("please specify big.matrix.file.path")

        }

        # convert to big.matrix
        dimnames(father.genetic.data) <- NULL
        father.bm <- as.big.matrix(father.genetic.data, type = "short", backingfile = "father_bm",
                                   backingpath = big.matrix.file.path, descriptorfile = "father_bm_desc.rds",
                                   binarydescriptor = TRUE)
        rm(father.genetic.data)

    } else if (!is.null(father.genetic.data) & class(father.genetic.data) == "big.matrix"){

        if (! describe(father.genetic.data)@description$type %in% c("short")){

            stop("father.genetic.data must be a big.matrix of type short. To convert, see function deepcopy from package bigmemory.")

        }

        if (describe(father.genetic.data)@description$sharedType != "FileBacked"){

            stop("father.genetic.data must be a file backed big.matrix (father.genetic.data@description$sharedType == 'FileBacked')")

        }

        father.bm <- father.genetic.data

    }


    # further split the input data by exposure, if specified
    if (!is.null(complement.genetic.data)){

        bm.list <- list(case = case.bm, complement = comp.bm)
        overall.bm.desc.list <- list(case = describe(case.bm),
                                     complement = describe(comp.bm))

    } else {

        bm.list <- list(case = case.bm, mother = mother.bm,
                        father = father.bm)
        overall.bm.desc.list <- list(case = describe(case.bm),
                                     mother = describe(mother.bm),
                                     father = describe(father.bm))

    }

    if (is.null(exposure)){

        exposure.levels <- NULL
        exposure.risk.levels <- NULL
        bm.desc.list <- overall.bm.desc.list


    } else {

        if (is.null(exposure.risk.levels)){

            exposure.risk.levels <- rep(1, length(unique(exposure)))

        } else {

            exposure.risk.levels <- unlist(exposure.risk.levels[as.character(unique(exposure))])

        }

        exposure.levels <- unique(exposure)
        storage.mode(exposure.levels) <- "integer"
        storage.mode(exposure.risk.levels) <- "integer"

        # list of rows for each exposure
        splits <- lapply(exposure.levels, function(exposure.level){

            exposure == exposure.level

        })

        # make big.matrix for each exposure level
        bm.desc.list <- lapply(seq_along(splits), function(split.number){

            these.rows <- splits[[split.number]]
            if (!is.null(complement.genetic.data)){

                case.bm.exp <- deepcopy(case.bm, rows = these.rows, type = "short",
                                      backingfile = paste0("case_bm_exp_", split.number),
                                      backingpath = big.matrix.file.path,
                                      descriptorfile = paste0("case_bm_exp_", split.number, "_desc.rds"),
                                      binarydescriptor = TRUE)

                comp.bm.exp <- deepcopy(comp.bm, rows = these.rows, type = "short",
                                        backingfile = paste0("comp_bm_exp_", split.number),
                                        backingpath = big.matrix.file.path,
                                        descriptorfile = paste0("comp_bm_exp_", split.number, "_desc.rds"),
                                        binarydescriptor = TRUE)

                bm.desc.list.exp <- list(case = describe(case.bm.exp),
                                         complement = describe(comp.bm.exp))

            } else {

                case.bm.exp <- deepcopy(case.bm, rows = these.rows, type = "short",
                                        backingfile = paste0("case_bm_exp_", split.number),
                                        backingpath = big.matrix.file.path,
                                        descriptorfile = paste0("case_bm_exp_", split.number, "_desc.rds"),
                                        binarydescriptor = TRUE)

                mother.bm.exp <- deepcopy(mother.bm, rows = these.rows, type = "short",
                                        backingfile = paste0("mother_bm_exp_", split.number),
                                        backingpath = big.matrix.file.path,
                                        descriptorfile = paste0("mother_bm_exp_", split.number, "_desc.rds"),
                                        binarydescriptor = TRUE)

                father.bm.exp <- deepcopy(father.bm, rows = these.rows, type = "short",
                                          backingfile = paste0("father_bm_exp_", split.number),
                                          backingpath = big.matrix.file.path,
                                          descriptorfile = paste0("father_bm_exp_", split.number, "_desc.rds"),
                                          binarydescriptor = TRUE)

                bm.desc.list.exp <- list(case = describe(case.bm.exp),
                                         mother = describe(mother.bm.exp),
                                         father = describe(father.bm.exp))

            }

            return(bm.desc.list.exp)

        })
        n.exposures <- length(bm.desc.list)
        bm.desc.list[[n.exposures + 1]] <- overall.bm.desc.list

    }

    # if needed compute sampling probs
    n.candidate.snps <- ncol(case.bm)

    if (is.null(snp.sampling.probs)){

        ### use conditional logistic regression to estimate univariate association ###

        if (is.null(categorical.exposures)){

            res.list <- bplapply(seq_len(n.candidate.snps), function(snp, bm.list) {

                case.snp <- bm.list$case[ , snp]
                missing.geno <- case.snp == -9
                case.snp <- case.snp[!missing.geno]
                n.fam <- length(case.snp)
                case.status <- c(rep(1, n.fam), rep(0, n.fam))
                ids <- rep(seq_len(n.fam), 2)

                if (length(bm.list) == 3){

                    mom.snp <- bm.list$mother[ , snp]
                    mom.snp <- mom.snp[!missing.geno]
                    dad.snp <- bm.list$father[ , snp]
                    dad.snp <- dad.snp[!missing.geno]
                    comp.snp <- mom.snp + dad.snp - case.snp

                } else {

                    comp.snp <- bm.list$complement[ , snp]
                    comp.snp <- comp.snp[!missing.geno]

                }

                # get p-value of association from conditional logistic regression
                case.comp.geno <- c(case.snp, comp.snp)
                clogit.res <- clogit(case.status ~ case.comp.geno + strata(ids), method = "approximate")
                clogit.chisq <- summary(clogit.res)$logtest[1]

                return(list(case.snp = case.snp, comp.snp = comp.snp, chisq = clogit.chisq))

            }, bm.list = bm.list, BPPARAM = bp.param)
            chisq.stats <- do.call("c", lapply(res.list, function(x) x$chisq))

        } else {

            res.list <- bplapply(seq_len(n.candidate.snps), function(snp, bm.list) {

                case.snp <- bm.list$case[ , snp]
                missing.geno <- case.snp == -9
                case.snp <- case.snp[!missing.geno]
                n.fam <- length(case.snp)
                exposure.var <- factor(rep(exposure[!missing.geno], 2))
                case.status <- c(rep(1, n.fam), rep(0, n.fam))
                ids <- rep(seq_len(n.fam), 2)

                if (length(bm.list) == 3){

                    mom.snp <- bm.list$mother[ , snp]
                    mom.snp <- mom.snp[!missing.geno]
                    dad.snp <- bm.list$father[ , snp]
                    dad.snp <- dad.snp[!missing.geno]
                    comp.snp <- mom.snp + dad.snp - case.snp

                } else {

                    comp.snp <- bm.list$complement[ , snp]
                    comp.snp <- comp.snp[!missing.geno]

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

            }, bm.list = bm.list, BPPARAM = bp.param)
            chisq.stats <- do.call("c", lapply(res.list, function(x) x$chisq))

        }

    }

    # take cumulative sum of ld.block.vec for output
    out.ld.vec <- cumsum(ld.block.vec)
    storage.mode(out.ld.vec) <- "integer"

    #### clean up chisq stats for models that did not converge ###
    chisq.stats[chisq.stats <= 0] <- 10^-10
    chisq.stats[is.infinite(chisq.stats)] <- max(chisq.stats[is.finite(chisq.stats)])

    return(list(genetic.data.list = bm.desc.list, chisq.stats = chisq.stats, ld.block.vec = out.ld.vec,
        exposure = exposure, exposure.levels = exposure.levels, exposure.risk.levels = exposure.risk.levels))

}
