#' A function to pre-process case-parent triad or disease-discordant sibling data.
#'
#' This function performs several pre-processing steps, intended for use before
#' function run.gadgets.
#'
#' @param case.genetic.data The genetic data of the disease affected children
#' from case-parent trios or disease-discordant sibling pairs. If searching for
#' maternal SNPs that are related to risk of disease in the child, some of the
#' columns in \code{case.genetic.data} may contain maternal SNP genotypes
#' (See argument \code{mother.snps} for how to indicate which SNPs columns
#' correspond to maternal genotypes). Columns are SNP allele counts, and rows
#' are individuals. This object may either be of class matrix' OR of class
#' 'big.matrix'. If of class 'big.matrix' it must be file backed as type
#' 'integer' (see the \code{bigmemory} package for more information). The
#' ordering of the columns must be consistent with the LD structure specified
#' in \code{ld.block.vec}. The genotypes cannot be  dosages imputed with
#' uncertainty.
#' @param complement.genetic.data A genetic dataset for the controls
#' corresponding to the genotypes in \code{case.genetic.data}.For SNPs that
#' correspond to the affected child in \code{case.genetic.data}, the
#' corresponding column in \code{complement.genetic.data} should be set equal to
#' mother allele count + father allele count - case allele count. If using
#' disease-discordant siblings this argument should be the genotypes for the
#' unaffected siblings. For SNPs in \code{case.genetic.data} that represent
#' maternal genotypes (if any) the corresponding column in
#' \code{complement.genetic.data} should be the paternal genotypes for that SNP.
#' Regardless, \code{complement.genetic.data} may be an object of either class
#' matrix' OR of class 'big.matrix'. If of class 'big.matrix' it must be file
#' backed as type 'integer' (see the bigmemory package for more information).
#' Columns are SNP allele counts, rows are families. If not specified,
#' \code{father.genetic.data} and \code{mother.genetic.data} must be specified.
#' The genotypes cannot be dosages imputed with uncertainty.
#' @param father.genetic.data The genetic data for the fathers of the cases in
#' \code{case.genetic.data}. This should only be specified when searching for
#' epistasis or GxGxE effects based only on case-parent triads, and not when
#' searching for maternal SNPs that are related to the child's risk of disease.
#' Columns are SNP allele counts, rows are individuals. This object may either
#' be of class 'matrix' OR of class 'big.matrix'. If of class big.matrix' it
#' must be file backed as type 'integer' (see the bigmemory package for more
#' information). The genotypes cannot be dosages imputed with uncertainty.
#' @param mother.genetic.data The genetic data for the mothers of the cases in
#' \code{case.genetic.data}. This should only be specified when searching for
#' epistasis or GxGxE effects based only on case-parent triads, and not when
#' searching for maternal SNPs that are related to the child's risk of disease.
#' Columns are SNP allele counts, rows are individuals. This object may either
#' be of class 'matrix' OR of class 'big.matrix'. If of class big.matrix' it
#' must be file backed as type 'integer' (see the bigmemory package for more
#' information). The genotypes cannot be dosages imputed with uncertainty.
#' @param ld.block.vec An integer vector specifying the linkage blocks of the
#' input SNPs. As an example, for 100 candidate SNPs, suppose we specify
#' \code{ld.block.vec <- c(25, 50, 25)}. This vector indicates that the input
#' genetic data has 3 distinct linkage blocks, with SNPs 1-25 in the first
#' linkage block, 26-75 in the second block, and 76-100 in the third block.
#' Note that this means the ordering of the columns (SNPs) in
#' \code{case.genetic.data} must be consistent with the LD blocks specified in
#' \code{ld.block.vec}. In the absence of outside information, a reasonable
#' default is to consider SNPs to be in LD if they are located on the same
#' biological chromosome. If \code{case.genetic.data} includes both maternal
#' and child SNP genotypes, we recommend considering any maternal SNP and any
#' child SNP located on the same nominal biological chromosome as 'in linkage'.
#' E.g., we recommend considering any maternal SNPs located on chromosome 1
#' as being 'linked' to any child SNPs located on chromosome 1, even though,
#' strictly speaking, the maternal and child SNPs are located on separate pieces
#' of DNA. If not specified, \code{ld.block.vec} defaults to assuming all input
#' SNPs are in linkage, which may be overly conservative and could
#' adversely affect performance.
#' @param bp.param The BPPARAM argument to be passed to bplapply when
#' estimating marginal disease associations for each SNP. If using a cluster
#' computer, this parameter needs to be set with care. See
#' \code{BiocParallel::bplapply} for more details.
#' @param snp.sampling.probs A vector indicating the sampling probabilities of
#' the SNPs in \code{case.genetic.data}. SNPs will be sampled in the
#' genetic algorithm proportional to the values specified. If not specified, by
#' default, chi-square statistics of association will be computed for
#' each SNP, and sampling will be proportional to the square root of those
#' statistics. If user specified, the values of \code{snp.sampling.probs} need
#' not sum to 1, they just need to be positive real numbers. See argument
#' \code{prob} from function \code{sample} for more details.
#' @param categorical.exposures (experimental) A matrix or data.frame of
#' integers corresponding to categorical exposures corresponding to the cases in
#' \code{case.genetic.data}. Defaults to NULL, which will result in GADGETS
#' looking for epistatic interactions, rather than SNP by exposure interactions.
#' \code{categorical.exposures} should not be missing any data; families with
#' missing exposure data should be removed from the analysis prior to input.
#' @param continuous.exposures (experimental) A matrix or data.frame of numeric
#' values representing continuous exposures corresponding to the families in
#' \code{case.genetic.data}. Defaults to NULL, which will result in GADGETS
#' searching for epistatic interactions, rather than SNP by exposure
#' interactions.
#' \code{continuous.exposures} should not be missing any data; families with
#' missing exposure data should be removed from the analysis prior to input.
#' @param mother.snps If searching for maternal SNPs that are associated
#' with disease in the child, the indices of the maternal SNP columns in object
#' \code{case.genetic.data}. Otherwise does not need to be specified.
#' @param child.snps If searching for maternal SNPs that are associated
#' with disease in the child, the indices of the child SNP columns in object
#' \code{case.genetic.data}. Otherwise does not need to be specified.
#' @param lower.order.gxe (experimental) A boolean indicating whether, if
#' multiple exposures of interest are input, E-GADGETS should search for only
#' for genetic interactions with the joint combination of exposures
#' (i.e., GxGxExE interactions), or if it should additionally search for
#' lower-order interactions that involve subsets of the exposures that were
#' input (i.e., GxGxE in addition to GxGxExE).
#' The default, FALSE, restricts the search to GxGxExE interactions. Users
#' should be cautious about including large numbers of input exposures, and, if
#' they do, very cautious about setting this argument to TRUE.
#' @return A list containing the following:
#' \describe{
#'  \item{case.genetic.data}{A matrix of case/maternal genotypes.}
#'  \item{complement.genetic.data}{A matrix of complement/sibling/paternal
#'  genotypes. If running E-GADGETS, this is set to a 1x1 matrix whose 
#'  single entry is 0, and not used}
#'  \item{mother.genetic.data}{If running E-GADGETS, A matrix of maternal 
#'  genotypes, otherwise a 1x1 matrix whose 
#'  single entry is 0.0, and not used}
#'  \item{father.genetic.data}{If running E-GADGETS, A matrix of mpaternal 
#'  genotypes, otherwise a 1x1 matrix whose 
#'  single entry is 0.0, and not used}
#'  \item{chisq.stats}{A vector of chi-square statistics corresponding to
#'  marginal SNP-disease associations, if \code{snp.sampling.probs}
#'  is not specified, and \code{snp.sampling.probs} otherwise.}
#'  \item{ld.block.vec}{A vector eaul to \code{cumsum(ld.block.vec)}.}
#'  \item{exposure.mat}{A design matrix of the input categorical and continuous
#'  exposures, if specified. Otherwise NULL.}
#'  \item{E_GADGETS}{A boolean indicating whether a GxGxE search is desired.}
#'  \item{mother.snps}{A vector of the column indices of maternal SNPs in
#'  \code{case.genetic.data}, set to NULL if not applicable.}
#'  \item{child.snps}{A vector of the column indices of child SNPs in
#'  \code{case.genetic.data}, set to NULL if not applicable.}
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
#' res <- preprocess.genetic.data(case[, 1:10],
#'                                father.genetic.data = dad[ , 1:10],
#'                                mother.genetic.data = mom[ , 1:10],
#'                                ld.block.vec = c(10))
#'
#' @importFrom bigmemory as.big.matrix describe attach.big.matrix
#' @importFrom data.table data.table rbindlist setorder
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom survival clogit strata coxph Surv
#' @importFrom stats as.formula cov model.matrix
#' @export

preprocess.genetic.data <- function(case.genetic.data,
                                    complement.genetic.data = NULL,
                                    father.genetic.data = NULL,
                                    mother.genetic.data = NULL,
                                    ld.block.vec = NULL,
                                    bp.param = bpparam(),
                                    snp.sampling.probs = NULL,
                                    categorical.exposures = NULL,
                                    continuous.exposures = NULL,
                                    mother.snps = NULL,
                                    child.snps = NULL, lower.order.gxe = FALSE) {

    #make sure the ld.block.vec is correctly specified
    if (is.null(ld.block.vec)){

        ld.block.vec <- ncol(case.genetic.data)

    } else {

        if (sum(ld.block.vec) != ncol(case.genetic.data)){

            stop("sum(ld.block.vec) must be equal to ncol(case.genetic.data)")

        }

    }

    # make sure the appropriate genetic data is included
    if (is.null(complement.genetic.data) & (is.null(father.genetic.data) |
                                            is.null(mother.genetic.data))) {

        stop("Must include complement.genetic.data or both father.genetic.data and mother.genetic.data")

    }

    ### if environmental exposures are provided, check their formatting ###
    if (!is.null(categorical.exposures)){

        cat.exposure.df <- as.data.frame(categorical.exposures)
        cat.exposure.df <- data.frame(apply(cat.exposure.df, 2, as.factor),
                                      stringsAsFactors = TRUE)

        # identify families with missing exposure data
        missing.exposure <- rowSums(is.na(cat.exposure.df)) > 0
        if (sum(missing.exposure) > 0){

            stop(paste("Please remove", sum(missing.exposure),
                       "families from analysis due to missing case categorical exposure(s)"))

        }

    }

    if (!is.null(continuous.exposures)){

        cont.exposure.df <- as.data.frame(continuous.exposures)

        # identify families with missing exposure data
        missing.exposure <- rowSums(is.na(cont.exposure.df)) > 0
        if (sum(missing.exposure) > 0){

            stop(paste("Please remove", sum(missing.exposure), "families from analysis due to case missing continuous exposure(s)"))

        }


    }

    # make full exposure df
    if (!is.null(continuous.exposures) & !is.null(categorical.exposures)){

        exposure.df <- cbind(cat.exposure.df, cont.exposure.df)

    } else if (!is.null(continuous.exposures) & is.null(categorical.exposures)){

        exposure.df <- cont.exposure.df

    } else if (is.null(continuous.exposures) & !is.null(categorical.exposures)){

        exposure.df <- cat.exposure.df

    } else {

        exposure.df <- NULL

    }
    
    # make sure we have trios for E-GADGETS
    E_GADGETS <- FALSE
    if (!is.null(exposure.df)){
        
        E_GADGETS <- TRUE
        
        if (is.null(mother.genetic.data) | is.null(father.genetic.data)){
            
            stop("E-GADGETS requires triads")
            
        }
        
    }

    # check formatting of input data and, if necessary, create memory  mapped
    # files
    if (!inherits(case.genetic.data, "matrix") &
        !inherits(case.genetic.data, "big.matrix")){

        stop("case.genetic.data must be of class matrix or big.matrix")

    }

    if (inherits(case.genetic.data, "matrix")){

        if (!all(round(case.genetic.data) == case.genetic.data, na.rm = TRUE)){

            stop("case.genetic.data genotypes must be integers, not dosages imputed with uncertainty")

        }

        storage.mode(case.genetic.data) <- "integer"

        # convert to big.matrix
        dimnames(case.genetic.data) <- NULL
        case.bm <- as.big.matrix(case.genetic.data, type = "integer",
                                 shared = FALSE)
        rm(case.genetic.data)
        gc(verbose = FALSE)

    } else if (inherits(case.genetic.data, "big.matrix")){

        if (! describe(case.genetic.data)@description$type %in% c("integer")){

            stop("case.genetic.data must be a big.matrix of type integer. To convert, see function deepcopy from package bigmemory.")

        }

        if (describe(case.genetic.data)@description$sharedType != "FileBacked"){

            stop("case.genetic.data must be a file backed big.matrix (case.genetic.data@description$sharedType == 'FileBacked')")

        }

        case.bm <- case.genetic.data

    }

    if (!is.null(complement.genetic.data) &
        !inherits(complement.genetic.data, "matrix") &
        !inherits(complement.genetic.data, "big.matrix")){

        stop("complement.genetic.data must be of class matrix or big.matrix")

    }

    if (!is.null(complement.genetic.data) &
        inherits(complement.genetic.data, "matrix")){

        if (!all(round(complement.genetic.data) == complement.genetic.data,
                 na.rm = TRUE)){

            stop("complement.genetic.data genotypes must be integers,
                 not dosages imputed with uncertainty")

        }

        storage.mode(complement.genetic.data) <- "integer"

        # convert to big.matrix
        dimnames(complement.genetic.data) <- NULL
        comp.bm <- as.big.matrix(complement.genetic.data, type = "integer",
                                 shared = FALSE)
        rm(complement.genetic.data)
        gc(verbose = FALSE)

    } else if (!is.null(complement.genetic.data) &
               inherits(complement.genetic.data, "big.matrix")){

        if (! describe(complement.genetic.data)@description$type %in%
            c("integer")){

            stop("complement.genetic.data must be a big.matrix of type integer. To convert, see function deepcopy from package bigmemory.")

        }

        if (describe(complement.genetic.data)@description$sharedType !=
            "FileBacked"){

            stop("complement.genetic.data must be a file backed big.matrix (complement.genetic.data@description$sharedType == 'FileBacked')")

        }

        comp.bm <- complement.genetic.data

    }

    if (!is.null(mother.genetic.data) &
        !inherits(mother.genetic.data,  "matrix") &
        !inherits(mother.genetic.data, "big.matrix")){

        stop("mother.genetic.data must be of class matrix or big.matrix")

    }

    if (!is.null(mother.genetic.data) &
        inherits(mother.genetic.data, "matrix")){

        if (!all(round(mother.genetic.data) == mother.genetic.data,
                 na.rm = TRUE)){

            stop("mother.genetic.data genotypes must be integers, not dosages imputed with uncertainty")

        }

        storage.mode(mother.genetic.data) <- "integer"

        # convert to big.matrix
        dimnames(mother.genetic.data) <- NULL
        mother.bm <- as.big.matrix(mother.genetic.data,
                                   type = "integer", shared = FALSE)
        rm(mother.genetic.data)
        gc(verbose = FALSE)

    } else if (!is.null(mother.genetic.data) &
               inherits(mother.genetic.data, "big.matrix")){

        if (! describe(mother.genetic.data)@description$type %in% c("integer")){

            stop("mother.genetic.data must be a big.matrix of type integer. To convert, see function deepcopy from package bigmemory.")

        }

        if (describe(mother.genetic.data)@description$sharedType !=
            "FileBacked"){

            stop("mother.genetic.data must be a file backed big.matrix (mother.genetic.data@description$sharedType == 'FileBacked')")

        }

        mother.bm <- mother.genetic.data

    }

    if (!is.null(father.genetic.data) &
        !inherits(father.genetic.data, "matrix") &
        !inherits(father.genetic.data, "big.matrix")){

        stop("father.genetic.data must be of class matrix or big.matrix")

    }

    if (!is.null(father.genetic.data) &
        inherits(father.genetic.data, "matrix")){

        if (!all(round(father.genetic.data) ==
                 father.genetic.data, na.rm = TRUE)){

            stop("father.genetic.data genotypes must be integers, not dosages imputed with uncertainty")

        }
        storage.mode(father.genetic.data) <- "integer"

        # convert to big.matrix
        dimnames(father.genetic.data) <- NULL
        father.bm <- as.big.matrix(father.genetic.data,
                                   type = "integer", shared = FALSE)
        rm(father.genetic.data)
        gc(verbose = FALSE)

    } else if (!is.null(father.genetic.data) &
               inherits(father.genetic.data, "big.matrix")){

        if (! describe(father.genetic.data)@description$type %in% c("integer")){

            stop("father.genetic.data must be a big.matrix of type integer. To convert, see function deepcopy from package bigmemory.")

        }

        if (describe(father.genetic.data)@description$sharedType
            != "FileBacked"){

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

        ### use conditional logistic regression to estimate association ###
        n.fam <- nrow(case.bm)
        n.candidate.snps <- ncol(case.bm)
        case.status <- c(rep(1, n.fam), rep(0, n.fam))
        ids <- rep(seq_len(n.fam), 2)

        if (!E_GADGETS){

            res.list <- bplapply(seq_len(n.candidate.snps),
                                 function(snp, bm.list) {

                case.snp <- bm.list$case[ , snp]
                if (length(bm.list) == 3){

                    mom.snp <- bm.list$mother[ , snp]
                    dad.snp <- bm.list$father[ , snp]
                    comp.snp <- mom.snp + dad.snp - case.snp

                } else {

                    comp.snp <- bm.list$complement[ , snp]

                }

                # get test stat from conditional logistic regression
                case.comp.geno <- c(case.snp, comp.snp)
                clogit.res <- clogit(case.status ~ case.comp.geno + strata(ids),
                                     method = "approximate")
                clogit.chisq <- summary(clogit.res)$logtest[1]

                return(list(case.snp = case.snp, comp.snp = comp.snp,
                            chisq = clogit.chisq))

            }, bm.list = bm.list, BPPARAM = bp.param)
            chisq.stats <- do.call("c", lapply(res.list, function(x) x$chisq))

        } else {

            exposures <- rbind(exposure.df, exposure.df)

            res.list <- bplapply(seq_len(n.candidate.snps),
                                 function(snp, bm.list, exposures) {

                case.snp <- bm.list$case[ , snp]
                if (length(bm.list) == 3){

                    mom.snp <- bm.list$mother[ , snp]
                    dad.snp <- bm.list$father[ , snp]
                    comp.snp <- mom.snp + dad.snp - case.snp

                } else {

                    comp.snp <- bm.list$complement[ , snp]

                }

                # get test stat from conditional logistic regression
                case.comp.geno <- c(case.snp, comp.snp)
                geno.df <- data.frame(case.status = case.status,
                                      case.comp.geno = case.comp.geno,
                                      ids = ids)
                df <- cbind(geno.df, exposures)

                #make model formula
                exposure.vars <- colnames(exposures)
                if (!lower.order.gxe){

                  exposure.part <- paste0("case.comp.geno:",
                                          paste(exposure.vars, collapse = ":"))

                } else {

                  exposure.part <- paste0("case.comp.geno*",
                                          paste(exposure.vars, collapse = "*"))

                }

                model.this <- as.formula(
                    paste0("case.status ~ case.comp.geno + ",
                           exposure.part, " + strata(ids)"))
                full.model <- clogit(model.this, method = "approximate",
                                     data  = df)
                full.model.ll <- full.model$loglik[2]
                reduced.model <- clogit(case.status ~ case.comp.geno +
                                            strata(ids),
                                        method = "approximate", data  = df)
                reduced.model.ll <- reduced.model$loglik[2]
                clogit.chisq <- 2*(full.model.ll - reduced.model.ll)

                return(list(case.snp = case.snp, comp.snp = comp.snp,
                            chisq = clogit.chisq))

            }, bm.list = bm.list, exposures = exposures, BPPARAM = bp.param)
            chisq.stats <- do.call("c", lapply(res.list, function(x) x$chisq))

        }

    }

    # take cumulative sum of ld.block.vec for output
    out.ld.vec <- cumsum(ld.block.vec)
    storage.mode(out.ld.vec) <- "integer"

    #### clean up chisq stats for models that did not converge ###
    chisq.stats[chisq.stats <= 0] <- 10^-10
    inf.idx <- is.infinite(chisq.stats)
    finite.idx <- is.finite(chisq.stats)
    chisq.stats[inf.idx] <- max(chisq.stats[finite.idx])
    
    case.data <- case.bm[]
    if (any(! case.data %in% c(NA, 0, 1, 2))){
        
        stop("Miscoded case genotypes, genotypes must be coded NA, 0, 1, or 2")
        
    }

    if (!E_GADGETS){
        
        if (!"complement" %in% names(bm.list)){
            
            comp.data <- mother.bm[] + father.bm[] - case.bm[]
            
        } else {
            
            comp.data <- comp.bm[]
            
        }
        
        # confirm no miscoded genotypes
        if (any(! comp.data %in% c(NA, 0, 1, 2))){
            
            stop("Miscoded genotypes, genotypes must be coded NA, 0, 1, or 2")
            
        }
        
        # set missing to -9
        if (any(is.na(case.data)) | any(is.na(comp.data))){
            
            case.data[is.na(case.data) | is.na(comp.data)] <- -9
            comp.data[is.na(case.data) | is.na(comp.data)] <- -9
            
        }
        
        #data not used, but required as input to run.gadgets
        mother.data <- matrix(0.0, 1, 1)
        father.data <- matrix(0.0, 1, 1)
        exposure.mat <- matrix(0.0, 1, 1)
        
    } else {
        
        comp.data <- matrix(0, 1, 1)
        storage.mode(comp.data) <- "integer"
        mother.data <- mother.bm[] 
        father.data <- father.bm[]
        
        if (any(! mother.data %in% c(NA, 0, 1, 2))){
            
            stop("Miscoded mother genotypes, genotypes must be coded NA, 0, 1, or 2")
            
        }
        
        if (any(! father.data %in% c(NA, 0, 1, 2))){
            
            stop("Miscoded father genotypes, genotypes must be coded NA, 0, 1, or 2")
            
        }
        
        #handle missing genotypes
        missing.geno <- is.na(case.data) | is.na(mother.data) | 
            is.na(father.data)
        
        if (any(missing.geno)){
            
            case.data[missing.geno] <- -9
            mother.data[missing.geno] <- -9
            father.data[missing.geno] <- -9
            
        }
        
        # format exposure matrix 
        exposure.vars <- colnames(exposure.df)
        if (!lower.order.gxe){
            
            model.terms <- paste0("~", paste(exposure.vars, collapse = ":"))
            
        } else {
            
            model.terms <- paste0("~", paste(exposure.vars, collapse = "*"))
            
        }
        exposure.mat <- model.matrix(as.formula(model.terms) , data = exposure.df)
        attributes(exposure.mat)$assign <- NULL
        attributes(exposure.mat)$contrasts <- NULL
        exposure.mat <- exposure.mat[ , -1, drop = FALSE]
        exposure.mat <- exposure.mat + 0.0
        
    }
    
    return(list(case.genetic.data = case.data,
                complement.genetic.data = comp.data,
                mother.genetic.data = mother.data,
                father.genetic.data = father.data,
                chisq.stats = chisq.stats,
                ld.block.vec = out.ld.vec,
                exposure.mat = exposure.mat,
                E_GADGETS = E_GADGETS,
                mother.snps = mother.snps,
                child.snps = child.snps))

}
