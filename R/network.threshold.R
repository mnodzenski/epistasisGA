#' A function to identify a threshold fitness score above which chromosomes should be included in a network plot
#'
#' This function identifies a threshold fitness score above which chromosomes should be included in a network plot
#'
#' @param observed.results A data frame of results of GA runs for the observed data, containing only unique chromosomes. The first n columns are snp ids, the next n columns are there associated mean vector components.
#' @param permutation.list A list of data frames, each containing unique results for GA runs, for permutations of the observed data.
#' @param threshold.val The desired theshold. This function will return a fitness score such that \code{threshold.val}*100\%  of permutations have at least one chromosome exceeding this score.

#' @return A numeric fitness score such that \code{threshold.val}*100\%  of permutations have at least one chromosome exceeding this score.
#'
#' @importFrom  data.table rbindlist
#' @export

network.threshold <- function(observed.results, permutation.list, threshold.val){

  #error checking
  if (any(duplicated(observed.results$fitness.score))){

    stop("Observed results must included only unique chromosomes, duplicates detected.")

  }

  #observed fitness scores
  obs.scores <- observed.results$fitness.score
  obs.scores <- sort(obs.scores, decreasing = TRUE)
  n.obs.scores <- length(obs.scores)

  #permutation fitness scores
  perm.scores <- unique(rbindlist(permutation.list)$fitness.score)
  n.perm.scores <- length(perm.scores)

  #find the observed fitness score closest to producing threshold.val
  prop.higher <- rep(NA, n.obs.scores)
  for (i in 1:n.obs.scores){

    val <- obs.scores[i]
    prop.higher[i] <- sum(perm.scores > val)/n.perm.scores

    if (prop.higher <= threshold.val){

      fitness.score.threshold <- val
      prop.permutations.higher <- prop.higher[i - 1]
      break

    }

  }

  #return the fitness score threshold and proportion of permutations with at least one
  #fitness score higher than this value
  return(list(fitness.score.threshold = fitness.score.threshold, prop.permutations.higher = prop.permutations.higher))

}


