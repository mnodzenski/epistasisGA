#' A function to identify a threshold fitness score above which chromosomes should be included in a network plot
#'
#' This function identifies a threshold fitness score above which chromosomes should be included in a network plot
#'
#' @param observed.results A data frame of results of GA runs for the observed data, containing only unique chromosomes.
#' @param permutation.list A list of data frames, each containing unique results for GA runs, for permutations of the observed data.
#' @param threshold.val The desired theshold. This function will return a fitness score such that \code{threshold.val}*100\%  of permutations have at least one chromosome exceeding this score.

#' @return A numeric fitness score such that \code{threshold.val}*100\%  of permutation chromosomes exceed this score.
#'
#' data(case)
#' data(dad)
#' data(mom)
#' library(Matrix)
#' chrom.mat <- as.matrix(bdiag(list(matrix(rep(TRUE, 2500^2), nrow = 2500),
#'                               matrix(rep(TRUE, 2500^2), nrow = 2500),
#'                               matrix(rep(TRUE, 2500^2), nrow = 2500),
#'                               matrix(rep(TRUE, 2500^2), nrow = 2500))))
#'
#  ## preprocess data
#' pp.list <- preprocess.genetic.data(case[, 1:10], father.genetic.data = dad[ , 1:10],
#'                                mother.genetic.data = mom[ , 1:10],
#'                                chrom.mat = chrom.mat[ , 1:10])
#' ## run GA for observed data
#'
#' run.ga(pp.list, n.chromosomes = 5, chromosome.size = 3, results.dir = "tmp_2",
#'        cluster.type = "interactive", registryargs = list(file.dir = "tmp_reg", seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  combined.res3 <- combine.islands("tmp_2")
#'  unlink("tmp_reg", recursive = TRUE)
#'
#' ## create three permuted datasets
#' set.seed(1400)
#' perm.data.list <- permute.dataset(case[ , 1:10],
#'                                   father.genetic.data = dad[ , 1:10],
#'                                   mother.genetic.data = mom[ , 1:10],
#'                                   n.permutations = 3)
#'
#' ## pre-process permuted data
#' p1.list <- preprocess.genetic.data(perm.data.list[["permutation1"]]$case,
#'                                    complement.genetic.data = perm.data.list[["permutation1"]]$comp,
#'                                    chrom.mat = chrom.mat[ , 1:10])
#'
#' p2.list <- preprocess.genetic.data(perm.data.list[["permutation2"]]$case,
#'                                    complement.genetic.data = perm.data.list[["permutation2"]]$comp,
#'                                    chrom.mat = chrom.mat[ , 1:10])
#'
#' p3.list <- preprocess.genetic.data(perm.data.list[["permutation3"]]$case,
#'                                    complement.genetic.data = perm.data.list[["permutation3"]]$comp,
#'                                    chrom.mat = chrom.mat[ , 1:10])
#'
#' ## run GA for permuted data
#'
#' run.ga(p1.list, n.chromosomes = 5, chromosome.size = 3, results.dir = "p1_tmp_3",
#'        cluster.type = "interactive", registryargs = list(file.dir = "tmp_reg", seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  p1.combined.res3 <- combine.islands("p1_tmp_3")
#'  unlink("tmp_reg", recursive = TRUE)
#'
#' run.ga(p2.list, n.chromosomes = 5, chromosome.size = 3, results.dir = "p2_tmp_3",
#'        cluster.type = "interactive", registryargs = list(file.dir = "tmp_reg", seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  p2.combined.res3 <- combine.islands("p2_tmp_3")
#'  unlink("tmp_reg", recursive = TRUE)
#'
#' run.ga(p3.list, n.chromosomes = 5, chromosome.size = 3, results.dir = "p3_tmp_3",
#'        cluster.type = "interactive", registryargs = list(file.dir = "tmp_reg", seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  p3.combined.res3 <- combine.islands("p3_tmp_3")
#'  unlink("tmp_reg", recursive = TRUE)
#'
#'  ## get threshold fitness score
#'  network.threshold(combined.res3$unique.results,
#'                    list(p1.combined.res3$unique.results,
#'                         p2.combined.res3$unique.results,
#'                         p3.combined.res3$unique.results),
#'                    0.2)
#'
#' @importFrom  data.table rbindlist
#' @export

network.threshold <- function(observed.results, permutation.list, threshold.val){

  #error checking
  if (any(duplicated(observed.results$chromosome))){

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

    if (prop.higher[i] >= threshold.val){

      fitness.score.threshold <- val
      prop.permutations.higher <- prop.higher[i - 1]
      break

    }

  }

  #return the fitness score threshold and proportion of permutations with at least one
  #fitness score higher than this value
  return(list(fitness.score.threshold = fitness.score.threshold, prop.permutations.higher = prop.permutations.higher))

}


