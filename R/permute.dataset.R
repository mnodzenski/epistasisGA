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
#'                                ld.block.vec = c(10),
#'                                big.matrix.file.path = "tmp_bm")
#' set.seed(15)
#' perm.data.list <- permute.dataset(pp.list, "tmp_perm", n.permutations = 1)
#'
#' unlink(c('tmp_perm', 'tmp_bm'))
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
    genetic.data.list <- preprocessed.list$genetic.data.list

    ### permute the data ###
    n.families <- genetic.data.list[[1]]@description$totalRows
    if (is.null(preprocessed.list$exposure)){

        permuted.data.list <- bplapply(seq_len(n.permutations), function(permute, n.families, genetic.data.list) {

            # flip the case/complement status for these families
            flip.these <- seq_len(n.families)[as.logical(rbinom(n.families, 1, 0.5))]

            # make deep copies of the input data
            trio.study <- length(genetic.data.list) == 3
            perm.genetic.data <- lapply(genetic.data.list, function(in.data.desc){

                in.data <- attach.big.matrix(in.data.desc)
                if (trio.study){

                    # don't need a new copy of the parent data
                    if (grepl("mother|father", in.data.desc@description$filename)){

                        perm.data <- in.data

                    } else {

                        perm.data.name <- paste0(in.data.desc@description$filename, "_p", permute)
                        desc.data.name <- paste0(perm.data.name, "_desc.rds")
                        perm.data <- deepcopy(in.data, type = "double", backingfile = perm.data.name,
                                              backingpath = permutation.data.file.path,
                                              descriptorfile = desc.data.name,
                                              binarydescriptor = TRUE)
                    }

                }

                return(perm.data)

            })
            names(perm.genetic.data) <- names(genetic.data.list)
            perm.data.addresses <- lapply(perm.genetic.data, function(x) x@address)
            perm.data.desc <- lapply(perm.genetic.data, function(x) describe(x))
            names(perm.data.desc) <- names(genetic.data.list)
            names(perm.data.addresses) <- names(genetic.data.list)
            create_permuted_data(perm.data.addresses, flip.these, trio.study)
            return(perm.data.desc)

        }, n.families = n.families, genetic.data.list = genetic.data.list, BPPARAM = bp.param)

    } else {

        permuted.data.list <- lapply(seq_len(n.permutations), function(permute) {

            shuffled.order <- sample(seq_along(categorical.exposures), length(categorical.exposures))
            exposure.perm <- categorical.exposures[shuffled.order]
            out.file <- file.path(permutation.data.file.path, paste0("exposure.p", permute, ".rds"))
            saveRDS(exposure.perm, out.file)
            return(exposure.perm)

        })

    }
    return(permuted.data.list)

}
