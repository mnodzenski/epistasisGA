#' A function to run a genetic algorithm on a collection of candidate snp sets, to identify epistatic variants.
#'
#' This function runs a genetic algorithm on a collection of candidate snp sets, to identify epistatic variants.
#'
#' @param case.genetic.data A genetic dataset from cases (for a dichotomous trait). Columns are snps, and rows are individuals.
#' @param father.genetic.data The genetic data for the father of the case. Columns are snps, rows are individuals.
#' @param mother.genetic.data The genetic data for the mother of the case. Columns are snps, rows are individuals.
#' @param n.chromosomes A scalar indicating the number of candidate collections of snps to use in the GA.
#' @param chromosome.size The number of snps within each candidate solution.
#' @param dist.type A character string indicating the type of distance measurement. Type 'knn' performs k-nearest neighbors classifications, type 'paired' performs paired length classifications.
#' @param generations The maximum number of generations for which the GA will run. Defaults to 2000.
#' @param k A numeric scalar corresponding to the number of nearest neighbors required for computing the fitness score. Defaults to 10
#' @param correct.thresh A numeric scalar between 0 and 1 indicating the minimum proportion of of cases among the nearest neighbors for a given individual for that individual to be considered correctly classified. Defaults to 0.9.
#' @param gen.same.fitness The number of consecutive generations with the same fitness score required for algorithm termination.
#' @param tol The maximum absolute pairwise difference among the top fitness scores from the previous 500 generations considered to be sufficient to stop producing new generations.
#' @param n.top.chroms The number of top scoring chromosomes, according to fitness score, to return.
#' @param run.parallel Logical indicator of whether the fitness scores should be computed in parallel using bplapply from bioconductor. If FALSE, will use sapply instead. Defaults to FALSE.
#' @return A list, whose first element is a list of the top 100 scoring chromosomes, and second element is a vector of the corresponding fitness scores.
#'
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#'
#' ga.res <- ga(case, dad, mom, 7, 3, dist.type = 'knn', generations = 1, k = 10,
#'  correct.thresh = 0.9, tol = 10^-6, run.parallel = TRUE)
#'
#' @importFrom Rfast Dist
#' @importFrom BiocParallel bplapply
#' @export

ga <- function(case.genetic.data, father.genetic.data, mother.genetic.data,
               n.chromosomes, chromosome.size, dist.type,
               generations = 2000, k = 10, correct.thresh = 0.9, gen.same.fitness = 500, tol = 10^-6,
               n.top.chroms = 100, run.parallel = F){

  ### Compute the complement data ###
  complement.genetic.data <- father.genetic.data + mother.genetic.data - case.genetic.data

  ### Compute matrices of differences between cases and complements ###
  case.minus.comp <- case.genetic.data - complement.genetic.data
  case.comp.different <- case.minus.comp != 0

  ### Compute matrix of expected squared difference between case and complement given parent geno ###
  expected.squared.diffs <- ifelse(mother.genetic.data == 1 & father.genetic.data == 1, 2,
                                   ifelse(abs(mother.genetic.data - father.genetic.data) == 1, 1, 0))

  ### initialize groups of candidate solutions ###
  chromosome.list <- vector(mode = "list", length = n.chromosomes)
  all.snps.idx <- 1:ncol(case.genetic.data)
  for (i in 1:n.chromosomes){

    snp.idx <- sort(sample(all.snps.idx, chromosome.size, replace = F))
    all.snps.idx <- all.snps.idx[-snp.idx]
    chromosome.list[[i]] <- snp.idx

  }

  ### Begin iterating over generations ###
  generation <- 1
  fitness.score.mat <- matrix(rep(NA, generations*n.chromosomes), nrow = generations)
  top.fitness <- rep(0, generations)
  last.gens.equal <- F
  top.generation.chromosome <- vector(mode = "list", length = generations)
  chromosome.mat <- matrix(rep(NA, generations*n.chromosomes), nrow = generations)

  while (generation <= generations & !last.gens.equal){

    print(paste("generation", generation))
    ### 1. compute the fitness score for each set of candidate snps ###
    print("Step 1/9")

    if (run.parallel){

      fitness.scores <- unlist(bplapply(1:length(chromosome.list), function(x) {

        fitness.score(case.genetic.data, complement.genetic.data, case.comp.different,
                      chromosome.list[[x]], dist.type, case.minus.comp, expected.squared.diffs,
                      k, correct.thresh)

      }))

    } else {

      fitness.scores <- sapply(1:length(chromosome.list), function(x) {

        fitness.score(case.genetic.data, complement.genetic.data, case.comp.different,
                      chromosome.list[[x]], dist.type, case.minus.comp, expected.squared.diffs,
                      k, correct.thresh)

      })

    }

    #store the fitness scores and elements (snps) of the chromosomes
    fitness.score.mat[generation, ] <- fitness.scores
    chromosome.mat[generation, ] <- sapply(chromosome.list, function(x) paste(x, collapse = "."))

    ### 2. identify the top scoring candidate solution(s) and fitness score ###
    print("Step 2/9")
    max.fitness <- max(fitness.scores)
    top.chromosome.idx <- which(fitness.scores == max.fitness)
    top.chromosomes <- chromosome.list[top.chromosome.idx]

    #if we have the same chromosome multiple times, just take one of them
    #and put the duplicates in the pool to be resampled
    if (length(unique(top.chromosomes)) == 1){

      top.chromosome <- top.chromosomes[1]
      duplicate.top.chromosomes <- top.chromosome.idx[-1]

    #if we have everyone with the same score, just choose the first one and
    #re-sample all other chromosomes
    } else if (length(unique(top.chromosomes)) == n.chromosomes){

      top.chromosome <- top.chromosomes[1]
      duplicate.top.chromosomes <- top.chromosome.idx[-1]

    } else {

      top.chromosome <- top.chromosomes
      duplicate.top.chromosomes <- NULL

    }

    ### 3. identify the lower scoring chromosomes ###
    print("Step 3/9")
    lower.chromosomes <- c(which(fitness.scores != max.fitness), duplicate.top.chromosomes)

    ### 4. Sample with replacement from the lower scoring chroms ###
    print("Step 4/9")
    sampled.lower.idx <- sample(lower.chromosomes, length(lower.chromosomes),
                                replace = T, prob = 1 + fitness.scores[lower.chromosomes])
    sampled.lower.chromosomes <- chromosome.list[sampled.lower.idx]

    ### 5. Determine whether each lower chromosome will be subject to mutation or crossing over ###
    print("Step 5/9")

    # only allowing cross-overs between distinct chrosomes
    # (i.e, if a chromosome was sampled twice, it can't cross over with itself)
    unique.lower.idx <- unique(sampled.lower.idx)

    #note: need at least two crossovers assigned, and need an even number
    cross.overs <- rep(F, length(unique.lower.idx))
    while (sum(cross.overs) < 2 | sum(cross.overs) %% 2 != 0){

      cross.overs <- rbinom(length(unique.lower.idx), 1, 0.5) == 1

    }

    #those not getting crossover will be mutated
    mutations <- !cross.overs

    ### 6. Execute crossing over for the relevant chromosomes ###
    print("Step 6/9")
    cross.over.positions <- match(unique.lower.idx[cross.overs], sampled.lower.idx)
    cross.over.starts <- seq(1, length(cross.over.positions), by = 2)
    for (i in cross.over.starts){

      #grab pair of chromosomes to cross over
      chrom1 <- sampled.lower.chromosomes[[cross.over.positions[i]]]
      chrom2 <- sampled.lower.chromosomes[[cross.over.positions[i+1]]]

      #check for overlapping snps
      matching.snp.positions <- match(chrom1, chrom2)
      matching.snp.positions <-  matching.snp.positions[!is.na(matching.snp.positions)]
      not.matching.snp.positions <- setdiff(1:4, match(a,b))

      #order the second snp, first by the overlapping snps and then randomly afterwards
      chrom2 <- chrom2[c(matching.snp.positions, not.matching.snp.positions)]

      #randomly sample a crossing over point, making sure we do not allow duplicate snps
      #and also making sure we don't simply swap chromosomes
      possible.cut.points <- which(chrom1 != chrom2)
      if (length(possible.cut.points) == chromosome.size){
        possible.cut.points <- possible.cut.points[-1]
      }
      if (all(chrom1 == chrom2)){

        cutpoint <- 1

      } else {

        cut.point <- sample(possible.cut.points, 1)

      }

      #cross
      chrom1.cross <- chrom1
      chrom1.cross[cut.point:chromosome.size] <- chrom2[cut.point:chromosome.size]
      chrom2.cross <- chrom2
      chrom2.cross[cut.point:chromosome.size] <- chrom1[cut.point:chromosome.size]

      #replace in the chromosome list
      sampled.lower.chromosomes[[cross.over.positions[i]]] <- chrom1.cross
      sampled.lower.chromosomes[[cross.over.positions[i+1]]] <- chrom2.cross

    }

    ### 7. Mutate the chromosomes that were not crossed over ###
    print("Step 7/9")
    mutation.positions <- (1:length(sampled.lower.chromosomes))[-cross.over.positions]
    for (i in mutation.positions){

      #grab the chromosome
      target.chrom <- sampled.lower.chromosomes[[i]]

      #determine which snps to mutate
      mutate.these <- rbinom(length(target.chrom), 1, 0.5) == 1

      if (any(mutate.these)){

        #remove the chromosome's snps from the pool of available snps
        #and sample new snps for the mutations
        mutated.snps <- sample((1:ncol(case.genetic.data))[-target.chrom], sum(mutate.these))

        #substitute in mutations
        sampled.lower.chromosomes[[i]][mutate.these] <- mutated.snps

      }

    }

    ### 8. Combine into new population (i.e., the final collection of chromosomes for the next generation)
    print("Step 8/9")
    chromosome.list <- lapply(c(top.chromosome, sampled.lower.chromosomes), sort)

    ### 9.Increment Iterators
    print("Step 9/9")
    top.fitness[generation] <- max.fitness
    top.generation.chromosome[[generation]] <- top.chromosome[[1]]
    print(paste0("Max fitness score:", max.fitness))
    print("Top Chromosome(s):")
    print(top.chromosome)
    if (generation >= gen.same.fitness){

      last.gens <- top.fitness[(generation - (gen.same.fitness -1)):generation]
      last.gens.equal <- abs(max(last.gens) - min(last.gens)) < tol

    }
    generation <- generation + 1

  }
  ### Return the best chromosomes and fitness scores ###
  last.generation <- generation - 1
  chromosome.vec <- as.vector(chromosome.mat)
  duplicated.chroms <- duplicated(chromosome.vec)
  unique.chromosome.vec <- chromosome.vec[!duplicated.chroms]
  fitness.score.vec <- as.vector(fitness.score.mat)[!duplicated.chroms]
  ordered.fitness.scores <- order(fitness.score.vec, decreasing = T)

  top.fitness.scores <- fitness.score.vec[ordered.fitness.scores][1:n.top.chroms]
  top.chroms <- unique.chromosome.vec[ordered.fitness.scores][1:n.top.chroms]

  return(list(chromosomes = top.chroms, fitness.scores = top.fitness.scores))

}





