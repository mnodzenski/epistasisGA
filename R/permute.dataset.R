#' A function to create permuted datasets for permutation based hypothesis testing.
#'
#' This function creates permuted datasets for permutation based hypothesis testing of GADGETS fitness scores.
#'
#' @param case.genetic.data The genetic data of the disease affected children. Columns are SNP allele counts, and rows are individuals.
#' @param complement.genetic.data A genetic dataset from the complements of the cases, where
#' \code{complement.genetic.data} = mother SNP counts + father SNP counts - case SNP counts. If using affected/unaffected
#' sibling pairs, this should be the genetic data for the unaffected sibling.
#' Columns are SNP allele counts, rows are families. If not specified, \code{father.genetic.data} and \code{mother.genetic.data} must be specified.
#' @param father.genetic.data The genetic data for the fathers of the cases. Columns are SNP allele counts, rows are individuals.
#' Does not need to be specified if \code{complement.genetic.data} is specified.
#' @param mother.genetic.data The genetic data for the mothers of the cases. Columns are SNP allele counts, rows are individuals.
#' Does not need to be specified if \code{complement.genetic.data} is specified.
#' @param n.permutations The number of permuted datasets to create.
#' @return A list of \code{n.permutations} pairs of case and complement data, where the observed case/complement status has been randomly flipped or not flipped.
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
permute.dataset <- function(case.genetic.data, complement.genetic.data = NULL, father.genetic.data = NULL,
    mother.genetic.data = NULL, n.permutations = 100) {

    # make sure the appropriate genetic data is included
    if (is.null(complement.genetic.data) & is.null(father.genetic.data) & is.null(mother.genetic.data)) {

        stop("Must include complement.genetic.data or both father.genetic.data and mother.genetic.data")

    }

    ### Compute the complement data if not provided ###
    if (!is.null(father.genetic.data) & !is.null(mother.genetic.data)) {

        complement.genetic.data <- father.genetic.data + mother.genetic.data - case.genetic.data

    }

    ### permute the data ###
    n.families <- nrow(case.genetic.data)
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
    permuted.data.list

}
