#' A function to create permutation datasets for permutation based hypothesis testing.
#'
#' This function creates permutation datasets for permutation based hypothesis testing.
#'
#' @param case.genetic.data A genetic dataset from cases (for a dichotomous trait). Columns are snps, and rows are individuals.
#' @param complement.genetic.data A genetic dataset from the complements of the cases, where \code{complement.genetic.data} = mother snp counts + father snp counts - case snp counts. Columns are snps, rows are families. If not specified, \code{father.genetic.data} and \code{mother.genetic.data} must be specified.
#' @param father.genetic.data The genetic data for the father of the case. Columns are snps, rows are individuals. Does not need to be specified if \code{complement.genetic.data} is specified.
#' @param mother.genetic.data The genetic data for the mother of the case. Columns are snps, rows are individuals. Does not need to be specified if \code{complement.genetic.data} is specified.
#' @param n.permutations The number of permuted datasets to create.
#' @param seed.val The starting seed for the random permutations of the data.
#' @return A list of \code{n.permutations} pairs of case and complement data, where the case complement status has been randomly flipped or not flipped.
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' perm.data.list <- permute.dataset(case, father.genetic.data = dad, mother.genetic.data = mom,
#'                  n.permutations = 100)
permute.dataset <- function(case.genetic.data, complement.genetic.data = NULL, father.genetic.data = NULL,
                                    mother.genetic.data = NULL, n.permutations = 100, seed.val = 1){

  #make sure the appropriate genetic data is included
  if (is.null(complement.genetic.data) & is.null(father.genetic.data) & is.null(mother.genetic.data)){

    stop("Must include complement.genetic.data or both father.genetic.data and mother.genetic.data")

  }

  ### Compute the complement data if not provided ###
  if (!is.null(father.genetic.data) & !is.null(mother.genetic.data)){

    complement.genetic.data <- father.genetic.data + mother.genetic.data - case.genetic.data

  }

  ### permute the data ###
  set.seed(seed.val)
  n.families <- nrow(case.genetic.data)
  permuted.data.list <- lapply(1:n.permutations, function(x){

    perm <- as.logical(rbinom(n.families,1,0.5))
    case.perm <- case.genetic.data
    comp.perm <- complement.genetic.data
    comp.copy <- complement.genetic.data
    case.copy <- case.genetic.data
    case.perm[perm, ] <- comp.copy[perm, ]
    comp.perm[perm, ] <- case.copy[perm, ]
    list(case = case.perm, comp = comp.perm)

  })
  names(permuted.data.list) <- paste0("permutation", 1:n.permutations)
  permuted.data.list

}
