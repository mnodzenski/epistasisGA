#' A function to run a genetic algorithm on a collection of candidate snp sets, to identify epistatic variants.
#'
#' This function runs a genetic algorithm on a collection of candidate snp sets, to identify epistatic variants.
#'
#' @param case.genetic.data A genetic dataset from cases (for a dichotomous trait). Columns are snps, and rows are individuals.
#' @param complement.genetic.data A genetic dataset from the complements of the cases, where \code{complement.genetic.data} = mother snp counts + father snp counts - case snp counts. Columns are snps, rows are families. If not specified, \code{father.genetic.data} and \code{mother.genetic.data} must be specified.
#' @param father.genetic.data The genetic data for the father of the case. Columns are snps, rows are individuals. Does not need to be specified if \code{complement.genetic.data} is specified.
#' @param mother.genetic.data The genetic data for the mother of the case. Columns are snps, rows are individuals. Does not need to be specified if \code{complement.genetic.data} is specified.
#' @param n.chromosomes A scalar indicating the number of candidate collections of snps to use in the GA.
#' @param seed.val An integer indicating the seed to be used for the random samples.
#' @param chromosome.size The number of snps within each candidate solution.
#' @param n.different.snps.weight The number by which the number different snps between case and control is multiplied in computing the family weights. Defaults to 2.
#' @param n.both.one.weight The number by which the number of different snps equal to 1 in both case and control is multiplied in computing the family weights. Defaults to 1.
#' @param weight.function A function that takes the weighted sum of the number of different snps and snps both equal to one as an argument, and returns a family weight. Defaults to the identity function.
#' @param min.allele.freq The minimum minor allele frequency in the parents required for a snp to be considered as a potential GA solution. Any snps with MAF < \code{min.allele.freq} in the parents will be omitted. Defaults to 0.01.
#' @param generations The maximum number of generations for which the GA will run. Defaults to 2000.
#' @param gen.same.fitness The number of consecutive generations with the same fitness score required for algorithm termination.
#' @param min.n.risk.set A scalar indicating the minimum number of individuals whose case - control difference vector must have sign consistent with the sign of the weighted sum of the differences vectors across families. Defaults to 10.
#' @param tol The maximum absolute pairwise difference among the top fitness scores from the previous 500 generations considered to be sufficient to stop producing new generations.
#' @param n.top.chroms The number of top scoring chromosomes, according to fitness score, to return.
#' @return A list, whose first element is a data.table of the top \code{n.top.chroms scoring chromosomes}, their fitness scores, and their difference vectors. The second element is a scalar indicating the number of generations required to identify a solution, and the third element is the number of snps filtered due to MAF < \code{min.allele.freq}.
#'
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#'
#' ga.res <- run.ga(case, dad, mom, 7, 3, seed.val = 10, generations = 1)
#'
#' @importFrom matrixStats colSds
#' @importFrom data.table data.table rbindlist setorder
#' @export

run.ga <- function(case.genetic.data, complement.genetic.data = NULL, father.genetic.data = NULL, mother.genetic.data = NULL,
                   n.chromosomes, chromosome.size, seed.val,n.different.snps.weight = 2, n.both.one.weight = 1,
                   weight.function = identity, min.allele.freq = 0.01, generations = 2000, gen.same.fitness = 500,
                   min.n.risk.set = 10, tol = 10^-6, n.top.chroms = 100){

  #make sure the appropriate genetic data is included
  if (is.null(complement.genetic.data) & is.null(father.genetic.data) & is.null(mother.genetic.data)){

    stop("Must include complement.genetic.data or both father.genetic.data and mother.genetic.data")

  }

  #set seed for reproducibility
  set.seed(seed.val)
  print(paste("Starting GA. Seed value:", seed.val))

  ### find the snps with MAF < minimum threshold in the cases ###
  if (!is.null(father.genetic.data) & !is.null(mother.genetic.data)){

    alt.allele.freqs <- colSums(father.genetic.data + mother.genetic.data)/(4*nrow(father.genetic.data))
    below.maf.threshold <- alt.allele.freqs > (1 - min.allele.freq) | alt.allele.freqs < min.allele.freq
    original.col.numbers <- which(!below.maf.threshold)
    names(original.col.numbers) <- NULL

    ### remove the snps not meeting the required allele frequency threshold ###
    father.genetic.data <- father.genetic.data[ , !below.maf.threshold]
    mother.genetic.data <- mother.genetic.data[ , !below.maf.threshold]
    case.genetic.data <- case.genetic.data[ , !below.maf.threshold]

    ### Compute the complement data ###
    complement.genetic.data <- father.genetic.data + mother.genetic.data - case.genetic.data

  } else if (!is.null(complement.genetic.data)){

    alt.allele.freqs <- colSums(case.genetic.data + complement.genetic.data)/(4*nrow(case.genetic.data))
    below.maf.threshold <- alt.allele.freqs > (1 - min.allele.freq) | alt.allele.freqs < min.allele.freq
    original.col.numbers <- which(!below.maf.threshold)
    names(original.col.numbers) <- NULL

    ### remove the snps not meeting the required allele frequency threshold ###
    case.genetic.data <- case.genetic.data[ , !below.maf.threshold]
    complement.genetic.data <- complement.genetic.data[ , !below.maf.threshold]

  }

  ### Compute matrices of differences between cases and complements ###
  case.minus.comp <- as.matrix(case.genetic.data - complement.genetic.data)
  case.comp.different <- case.minus.comp != 0

  ### Compute matrix of snp 'Z-scores' ###
  mean.snp.diffs <- colMeans(case.minus.comp)
  sd.snp.diffs <- colSds(case.minus.comp)
  snp.zscores <- mean.snp.diffs/sd.snp.diffs

  ### Compute matrix indicating whether both the case and control have 1 copy of the alt allele ###
  both.one.mat <- complement.genetic.data == 1 & case.genetic.data == 1

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
  chromosome.mat.list <- vector(mode = "list", length = generations)
  sum.dif.vec.list <- vector(mode = "list", length = generations)

  while (generation <= generations & !last.gens.equal){

    print(paste("generation", generation))
    ### 1. compute the fitness score for each set of candidate snps ###
    print("Step 1/9")

    fitness.score.list <- lapply(1:length(chromosome.list), function(x) {

        chrom.fitness.score(case.comp.different, chromosome.list[[x]], case.minus.comp, both.one.mat,
                            n.different.snps.weight, n.both.one.weight, weight.function, min.n.risk.set)

    })

    fitness.scores <- sapply(fitness.score.list, function(x) x$fitness.score)
    sum.dif.vecs <- t(sapply(fitness.score.list, function(x) x$sum.dif.vecs))

    #store the fitness scores, elements (snps) of the chromosomes, sum of the difference vectors
    fitness.score.mat[generation, ] <- fitness.scores
    chromosome.mat.list[[generation]] <- data.table(t(sapply(chromosome.list, function(x) original.col.numbers[x])))
    sum.dif.vec.list[[generation]] <- data.table(sum.dif.vecs)

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

    ### 4. Sample with replacement from the existing chromosomes ###
    print("Step 4/9")
    #allow the top scoring chromosome to be sampled, but only sample from the unique chromosomes available
    sample.these <- !duplicated(chromosome.list)
    sampled.lower.idx <- sample(which(sample.these), length(lower.chromosomes),
                                replace = T, prob = fitness.scores[sample.these])
    sampled.lower.chromosomes <- chromosome.list[sampled.lower.idx]
    sampled.lower.dif.vecs <- sum.dif.vecs[sampled.lower.idx , ]
    sampled.lower.fitness.scores <- fitness.scores[sampled.lower.idx]

    ### 5. Determine whether each lower chromosome will be subject to mutation or crossing over ###
    print("Step 5/9")

    # only allowing cross-overs between distinct chromosomes
    # (i.e, if a chromosome was sampled twice, it can't cross over with itself)
    unique.lower.idx <- unique(sampled.lower.idx)

    cross.overs <- rep(F, length(unique.lower.idx))
    #note: need at least two crossovers assigned, and need an even number
    if (length(unique.lower.idx) > 1){

      while (sum(cross.overs) < 2 | sum(cross.overs) %% 2 != 0){

        cross.overs <- rbinom(length(unique.lower.idx), 1, 0.5) == 1

      }

    }

    #those not getting crossover will be mutated
    mutations <- !cross.overs

    ### 6. Execute crossing over for the relevant chromosomes ###
    print("Step 6/9")
    if (any(cross.overs)){

      cross.over.positions <- match(unique.lower.idx[cross.overs], sampled.lower.idx)
      cross.over.starts <- seq(1, length(cross.over.positions), by = 2)
      for (i in cross.over.starts){

        #grab pair of chromosomes to cross over
        chrom1 <- sampled.lower.chromosomes[[cross.over.positions[i]]]
        chrom1.dif.vecs <- sampled.lower.dif.vecs[cross.over.positions[i], ]
        chrom1.fitness.score <- sampled.lower.fitness.scores[cross.over.positions[i]]

        chrom2 <- sampled.lower.chromosomes[[cross.over.positions[i+1]]]
        chrom2.dif.vecs <- sampled.lower.dif.vecs[cross.over.positions[i+1], ]
        chrom2.fitness.score <- sampled.lower.fitness.scores[cross.over.positions[i+1]]

        if (any(duplicated(chrom1))){

          stop("Duplicate elements in chrom1")

        }

        if (any(duplicated(chrom2))){

          stop("Duplicate elements in chrom2")

        }

        #check for overlapping snps
        c1.c2.matching.snp.positions <- match(chrom1, chrom2)
        c1.c2.matching.snp.positions <-  c1.c2.matching.snp.positions[!is.na(c1.c2.matching.snp.positions)]
        c1.c2.not.matching.snp.positions <- setdiff(1:chromosome.size, c1.c2.matching.snp.positions)

        c2.c1.matching.snp.positions <- match(chrom2, chrom1)
        c2.c1.matching.snp.positions <-  c2.c1.matching.snp.positions[!is.na(c2.c1.matching.snp.positions)]
        c2.c1.not.matching.snp.positions <- setdiff(1:chromosome.size, c2.c1.matching.snp.positions)

        #for the non-matching snps, order by the decreasing magnitude of the difference vector in the chromsome with the higher fitness score
        #and order by the increasing magntidue of the difference vector in the chromosome with the lower fitness score
        #**ultimately will be used to substitute the higher magnitude elements in the lower scoring chromosome for the lower magnitude
        #elements in the higher scoring chromosome**
        if (chrom1.fitness.score >= chrom2.fitness.score){

          c1.c2.not.matching.snp.positions <- c1.c2.not.matching.snp.positions[order(abs(chrom2.dif.vecs[c1.c2.not.matching.snp.positions]))]
          c2.c1.not.matching.snp.positions <- c2.c1.not.matching.snp.positions[order(abs(chrom1.dif.vecs[c2.c1.not.matching.snp.positions]), decreasing = T)]

        } else{

          c1.c2.not.matching.snp.positions <- c1.c2.not.matching.snp.positions[order(abs(chrom2.dif.vecs[c1.c2.not.matching.snp.positions]), decreasing = T)]
          c2.c1.not.matching.snp.positions <- c2.c1.not.matching.snp.positions[order(abs(chrom1.dif.vecs[c2.c1.not.matching.snp.positions]))]

        }

        #order the chromosomes, first by the overlapping snps and the non-overlapping
        chrom2 <- chrom2[c(c1.c2.matching.snp.positions, c1.c2.not.matching.snp.positions)]
        chrom1 <- chrom1[c(c2.c1.matching.snp.positions, c2.c1.not.matching.snp.positions)]

        #determine how many snps could be crossed over
        #and also make sure we don't simply swap chromosomes
        possible.cut.points <- sort(which(chrom1 != chrom2), decreasing = T)
        if (length(possible.cut.points) == chromosome.size){

          n.possible.crosses <- chromosome.size - 1

        } else {

          n.possible.crosses <- length(possible.cut.points)

        }

        #determine how many snps will actually be crossed over
        n.crosses <- sample.int(n.possible.crosses, 1)

        #pick out their positions
        cross.points <- c(1:chromosome.size)[(chromosome.size - n.crosses + 1):chromosome.size]

        #exchange the high magnitude elements from the lower scoring chromosome with the low magnitude elements from the high
        #scoring chromsosome
        chrom1.cross <- chrom1
        chrom1.cross[cross.points] <- chrom2[cross.points]
        chrom2.cross <- chrom2
        chrom2.cross[cross.points] <- chrom1[cross.points]

        #error checking
        if(length(chrom1.cross) != chromosome.size){

          stop("chrom1 wrong length")
        }

        if(length(chrom2.cross) != chromosome.size){

          stop("chrom2 wrong length")
        }

        #replace in the chromosome list
        sampled.lower.chromosomes[[cross.over.positions[i]]] <- chrom1.cross
        sampled.lower.chromosomes[[cross.over.positions[i+1]]] <- chrom2.cross

      }

    }

    ### 7. Mutate the chromosomes that were not crossed over ###
    print("Step 7/9")
    mutation.positions <- (1:length(sampled.lower.chromosomes))[-cross.over.positions]
    for (i in mutation.positions){

      #grab the chromosome and its difference vector
      target.chrom <- sampled.lower.chromosomes[[i]]
      target.dif.vec <- as.vector(t(sampled.lower.dif.vecs[i, ]))

      #sort the chromosome elements from lowest absolute difference vector to highest
      target.chrom <- target.chrom[order(abs(target.dif.vec))]

      #determine which snps to mutate
      total.mutations <- sample.int(chromosome.size, 1)
      mutate.these <- 1:total.mutations

      #remove the chromosome's snps from the pool of available snps
      #and sample new snps for the mutations
      mutated.snps <- sample((1:ncol(case.genetic.data))[-target.chrom], total.mutations, prob = abs(snp.zscores)[-target.chrom])

      #substitute in mutations
      target.chrom[mutate.these] <- mutated.snps
      sampled.lower.chromosomes[[i]] <- target.chrom

    }

    ### 8. Combine into new population (i.e., the final collection of chromosomes for the next generation)
    print("Step 8/9")
    chromosome.list <- lapply(c(top.chromosome, sampled.lower.chromosomes), sort)

    ### 9.Increment Iterators
    print("Step 9/9")
    top.fitness[generation] <- max.fitness
    top.generation.chromosome[[generation]] <- original.col.numbers[top.chromosome[[1]]]
    print(paste0("Max fitness score:", max.fitness))
    print("Top Chromosome(s):")
    print(original.col.numbers[top.chromosome[[1]]])
    if (generation >= gen.same.fitness){

      last.gens <- top.fitness[(generation - (gen.same.fitness -1)):generation]
      last.gens.equal <- abs(max(last.gens) - min(last.gens)) < tol

    }
    generation <- generation + 1

  }
  ### Return the best chromosomes, their fitness scores, difference vectors and the number of generations ###
  last.generation <- generation - 1
  all.chrom.dt <- rbindlist(chromosome.mat.list[1:last.generation], use.names = F)
  all.chrom.dif.vec.dt <- rbindlist(sum.dif.vec.list[1:last.generation], use.names = F)
  unique.chromosome.dt <- unique(all.chrom.dt)
  colnames(unique.chromosome.dt) <- paste0("snp", 1:ncol(unique.chromosome.dt))
  unique.chrom.dif.vec.dt <- all.chrom.dif.vec.dt[!duplicated(all.chrom.dt), ]
  colnames(unique.chrom.dif.vec.dt) <- paste0("snp", 1:ncol(unique.chrom.dif.vec.dt), ".diff.vec")
  unique.fitness.score.vec <- as.vector(t(fitness.score.mat[1:last.generation, ]))[!duplicated(all.chrom.dt)]
  unique.results <- cbind(unique.chromosome.dt, unique.chrom.dif.vec.dt)
  unique.results[ , fitness.score := unique.fitness.score.vec]
  setorder(unique.results, -fitness.score)
  final.result <- unique.results[1:n.top.chroms, ]

  print(paste("Algorithm terminated after", last.generation, "runs."))
  return(list(top.chromosome.results = final.result, n.generations = last.generation, n.filtered.snps = sum(below.maf.threshold) ))

}


