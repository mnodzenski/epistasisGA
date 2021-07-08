#' A function to create permuted datasets for permutation based hypothesis testing.
#'
#' This function creates permuted datasets for permutation based hypothesis testing of GADGETS fitness scores.
#'
#' @param case.genetic.data The genetic data of the disease affected children. Columns are SNP allele counts, and rows are individuals. If running
#' permutations for a GxE search, this argument should not be specified.
#' @param complement.genetic.data A genetic dataset from the complements of the cases, where
#' \code{complement.genetic.data} = mother SNP counts + father SNP counts - case SNP counts. If using affected/unaffected
#' sibling pairs, this should be the genetic data for the unaffected sibling. If running
#' permutations for a GxE search, this argument should not be specified.
#' Columns are SNP allele counts, rows are families. If not specified, \code{father.genetic.data} and \code{mother.genetic.data} must be specified.
#' @param father.genetic.data The genetic data for the fathers of the cases. Columns are SNP allele counts, rows are individuals.
#' Does not need to be specified if \code{complement.genetic.data} is specified. If running
#' permutations for a GxE search, this argument should not be specified.
#' @param mother.genetic.data The genetic data for the mothers of the cases. Columns are SNP allele counts, rows are individuals.
#' Does not need to be specified if \code{complement.genetic.data} is specified. If running
#' permutations for a GxE search, this argument should not be specified.
#' @param categorical.exposures A vector of integers corresponding to categorical exposures for the cases. Defaults to NULL,
#' which will result in GADGETS looking for epistatic interactions, rather than SNP by exposure interactions. Any missing
#' exposure data should be coded as NA. Does not need to be specified unless a GxE search is being run.
#' @param n.permutations The number of permuted datasets to create.
#' @return If genetic data are specified, a list of \code{n.permutations} pairs of case and complement data,
#' where the observed case/complement status has been randomly flipped or not flipped. If exposure data are
#' specified a list of \code{n.permutations} vectors where the exposures have been randomly shuffled.
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' set.seed(15)
#' perm.data.list <- permute.dataset(case, father.genetic.data = dad, mother.genetic.data = mom,
#'                  n.permutations = 2)
#'
#' @export
permute.dataset <- function(case.genetic.data = NULL, complement.genetic.data = NULL, father.genetic.data = NULL,
    mother.genetic.data = NULL, categorical.exposures = NULL, n.permutations = 100) {

    # make sure the appropriate genetic data is included
    if (is.null(complement.genetic.data) & is.null(father.genetic.data) & is.null(mother.genetic.data)
        & !is.null(case.genetic.data)) {

        stop("Must include complement.genetic.data or both father.genetic.data and mother.genetic.data")

    }

    if (!is.null(categorical.exposures) & (!is.null(complement.genetic.data) | !is.null(father.genetic.data)
                                           | !is.null(mother.genetic.data) | !is.null(case.genetic.data))) {

        stop("If running permutations for GxE search, please do not specify any genetic.data arguments.")

    }

    ### Compute the complement data if not provided ###
    if (!is.null(case.genetic.data) & !is.null(father.genetic.data) & !is.null(mother.genetic.data)) {

        complement.genetic.data <- father.genetic.data + mother.genetic.data - case.genetic.data

    }

    ### permute the data ###
    n.families <- nrow(case.genetic.data)
    if (!is.null(case.genetic.data)){

        permuted.data.list <- lapply(seq_len(n.permutations), function(x) {

            perm <- as.logical(rbinom(n.families, 1, 0.5))
            case.perm <- case.genetic.data
            comp.perm <- complement.genetic.data
            comp.copy <- complement.genetic.data
            case.copy <- case.genetic.data
            case.perm[perm, ] <- comp.copy[perm, ]
            comp.perm[perm, ] <- case.copy[perm, ]
            list(case = case.perm, comp = comp.perm)

        })
        names(permuted.data.list) <- paste0("permutation", seq_len(n.permutations))

    } else {

        permuted.data.list <- lapply(seq_len(n.permutations), function(x) {

            shuffled.order <- sample(seq_along(categorical.exposures), length(categorical.exposures))
            exposure.perm <- categorical.exposures[shuffled.order]
            return(exposure.perm)


        })
        names(permuted.data.list) <- paste0("permutation", seq_len(n.permutations))

    }

    permuted.data.list

}
