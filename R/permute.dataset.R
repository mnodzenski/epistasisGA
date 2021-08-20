#' A function to create permuted datasets for permutation based hypothesis testing.
#'
#' This function creates permuted datasets for permutation based hypothesis testing of GADGETS fitness scores.
#'
#' @param preprocessed.list The output list from \code{preprocess.genetic.data} for the original genetic data.
#' @param permutation.matrix.file.path  If runing GADGETS for GxG interactions, this argument specifies a directory
#'  where memory mapped files of class 'big.memory' will be saved for each permuted dataset on disk. If searching
#'  for GxE interactions, permuted versions of the exposure vector will be saved to this directory.
#' @param n.permutations The number of permuted datasets to create.
#' @param bp.param The BPPARAM argument to be passed to bplapply when estimating marginal disease associations for each SNP.
#'  If using a cluster computer, this parameter needs to be set with care. See \code{BiocParallel::bplapply} for more details
#' @return If genetic data are specified, a list of \code{n.permutations} pairs of case and complement data,
#' where the observed case/complement status has been randomly flipped or not flipped. If exposure data are
#' specified a list of \code{n.permutations} vectors where the exposures have been randomly shuffled.
#' @examples
#'
#' data(case)
#' case <- as.matrix(case)
#' data(dad)
#' dad <- as.matrix(dad)
#' data(mom)
#' mom <- as.matrix(mom)
#' pp.list <- preprocess.genetic.data(case[, 1:10], father.genetic.data = dad[ , 1:10],
#'                                mother.genetic.data = mom[ , 1:10],
#'                                ld.block.vec = c(10))
#' set.seed(15)
#' perm.data.list <- permute.dataset(pp.list, "tmp_perm", n.permutations = 1)
#'
#' unlink(c('tmp_perm'))
#'
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom bigmemory deepcopy attach.big.matrix describe
#' @export
permute.dataset <- function(preprocessed.list, permutation.data.file.path, n.permutations = 100,
                            bp.param = bpparam()) {

    if (!dir.exists(permutation.data.file.path)){

        dir.create(permutation.data.file.path, recursive = TRUE)

    }
    permutation.data.file.path <- normalizePath(permutation.data.file.path)

    # grab input genetic data
    case.genetic.data <- preprocessed.list$case.genetic.data
    complement.genetic.data <- preprocessed.list$complement.genetic.data

    ### permute the data ###
    n.families <- nrow(case.genetic.data)
    if (is.null(preprocessed.list$exposure)){

        permuted.data.list <- bplapply(seq_len(n.permutations), function(permute, n.families, case.genetic.data,
                                                                         complement.genetic.data) {

            # flip the case/complement status for these families
            flip.these <- seq_len(n.families)[as.logical(rbinom(n.families, 1, 0.5))]
            case.perm <- case.genetic.data
            comp.perm <- complement.genetic.data
            case.perm[flip.these, ] <- complement.genetic.data[flip.these, ]
            comp.perm[flip.these, ] <- case.genetic.data[flip.these, ]
            case.out.file <- file.path(permutation.data.file.path, paste0("case.permute", permute, ".rds"))
            comp.out.file <- file.path(permutation.data.file.path, paste0("complement.permute", permute, ".rds"))
            saveRDS(case.perm, case.out.file)
            saveRDS(comp.perm, comp.out.file)

        }, n.families = n.families, case.genetic.data = case.genetic.data, complement.genetic.data = complement.genetic.data,
        BPPARAM = bp.param)

    } else {

        exposure <- preprocessed.list$exposure
        permuted.data.list <- lapply(seq_len(n.permutations), function(permute) {

            shuffled.order <- sample(seq_along(exposure), length(exposure))
            exposure.perm <- exposure[shuffled.order]
            out.file <- file.path(permutation.data.file.path, paste0("exposure.permute", permute, ".rds"))
            saveRDS(exposure.perm, out.file)

        })

    }

}
