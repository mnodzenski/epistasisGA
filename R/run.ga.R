#' A function to run a genetic algorithm on a collection of candidate snp sets, to identify epistatic variants.
#'
#' This function runs a genetic algorithm on a collection of candidate snp sets, to identify epistatic variants.
#'
#' @param data.list The output list from \code{preprocess.genetic.data}.
#' @param n.chromosomes A scalar indicating the number of candidate collections of snps to use in the GA.
#' @param chromosome.size The number of snps within each candidate solution.
#' @param results.dir The directory to which island results will be saved.
#' @param n.different.snps.weight The number by which the number different snps between case and control is multiplied in computing the family weights. Defaults to 2.
#' @param n.both.one.weight The number by which the number of different snps equal to 1 in both case and control is multiplied in computing the family weights. Defaults to 1.
#' @param weight.function A function that takes the weighted sum of the number of different snps and snps both equal to one as an argument, and returns a family weight. Defaults to the identity function.
#' @param generations The maximum number of generations for which the GA will run. Defaults to 2000.
#' @param gen.same.fitness The number of consecutive generations with the same fitness score required for algorithm termination.
#' @param tol The maximum absolute pairwise difference among the top fitness scores from the previous 500 generations considered to be sufficient to stop producing new generations.
#' @param n.top.chroms The number of top scoring chromosomes, according to fitness score, to return.
#' @param initial.sample.duplicates A logical indicating whether the same snp can appear in more than one chromosome in the initial sample of chromosomes (the same snp may appear in more than one chromosome thereafter, regardless). Default to F.
#' @param snp.sampling.type A string indicating how snps are to be sampled for mutations. Options are "zscore" or "random". Defaults to "zscore".
#' @param crossover.prop A numeric between 0 and 1 indicating the proportion of chromosomes to be subjected to cross over. The remaining proportion will be mutated. Defaults to 0.5.
#' @param n.islands An integer indicating the number of islands to be used.
#' @param island.cluster.size An integer equal to the number of islands among which population migration may occur.
#' @param migration.generations An integer equal to the number of generations between migration among islands.
#' @param n.migrations The number of chromosomes that migrate among islands.
#'
#' @return A list, whose first element is a data.table of the top \code{n.top.chroms scoring chromosomes}, their fitness scores, and their difference vectors. The second element is a scalar indicating the number of generations required to identify a solution.
#'
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' library(Matrix)
#' chrom.mat <- as.matrix(bdiag(list(matrix(rep(TRUE, 2500^2), nrow = 2500),
#'                               matrix(rep(TRUE, 2500^2), nrow = 2500),
#'                               matrix(rep(TRUE, 2500^2), nrow = 2500),
#'                               matrix(rep(TRUE, 2500^2), nrow = 2500))))
#' #ga.res <- run.ga(case, father.genetic.data = dad, mother.genetic.data = mom, n.chromosomes = 7,
#'  #                chromosome.size = 3, chrom.mat = chrom.mat, seed.val = 10, generations = 1)
#'
#' @importFrom matrixStats colSds rowMaxs
#' @importFrom data.table data.table rbindlist setorder
#' @importFrom stats rbinom sd
#' @importFrom survival clogit
#' @export

run.ga <- function(data.list, n.chromosomes, chromosome.size, results.dir,
                   n.different.snps.weight = 2, n.both.one.weight = 1,
                   weight.function = identity,generations = 500, gen.same.fitness = 50,
                   tol = 10^-6, n.top.chroms = 100, initial.sample.duplicates = F,
                   snp.sampling.type = "chisq", crossover.prop = 0.8, n.islands = 1000,
                   island.cluster.size = 4, migration.generations = 50, n.migrations = 20,
                   starting.seeds = NULL){

  ### make sure the island cluster size divides the number of islands evenly ###
  if ( n.islands %% island.cluster.size != 0){

    stop("Number of islands needs to be an multiple of island cluster size")

  }

  ### make sure, if starting seeds are provided, that enough are provided ###
  if (!is.null(starting.seeds)){

    if (length(starting.seeds) != n.islands){

      stop("Number of starting seeds is not equal to the number of islands")

    }

  }

  #### grab the analysis data ###
  case.genetic.data <- data.list$case.genetic.data
  complement.genetic.data <- data.list$complement.genetic.data
  original.col.numbers <- data.list$original.col.numbers
  chisq.stats <- data.list$chisq.stats
  chrom.mat <- data.list$chrom.mat

  #### clean up chisq stats for models that did not converge ###
  chisq.stats[chisq.stats <= 0] <- 10^-10
  chisq.stats[is.infinite(chisq.stats)] <- max(chisq.stats[is.finite(chisq.stats)])

  ### Compute matrices of differences between cases and complements ###
  case.minus.comp <- sign(as.matrix(case.genetic.data - complement.genetic.data))
  case.comp.different <- case.minus.comp != 0

  ### Compute matrix indicating whether both the case and control have 1 copy of the minor allele ###
  both.one.mat <- complement.genetic.data == 1 & case.genetic.data == 1

  ### set sampling type for mutation snps ###
  if (snp.sampling.type == "chisq"){

    snp.chisq <- sqrt(chisq.stats)

  } else if (snp.sampling.type == "random") {

    snp.chisq <- rep(1, ncol(case.minus.comp))

  }

  ### initialize seeds for different islands ###
  if (is.null(starting.seeds)){

    starting.seeds <- 1:n.islands

  }

  if (!dir.exists(results.dir)){

    dir.create(results.dir, recursive = T)

  } else {

    prev.islands <- list.files(results.dir)
    prev.islands <- prev.islands[! prev.islands %in% c("Old", "old") ]
    n.prev.islands <- length(prev.islands)
    n.islands <- n.islands - n.prev.islands
    prev.seeds <- as.numeric(gsub("island|.rds", "", prev.islands))
    starting.seeds <- setdiff(starting.seeds, prev.seeds)

    if (length(starting.seeds) == 0){

      stop("All islands have already been evolved")

    } else {

      starting.seeds <- setdiff(starting.seeds, prev.seeds)[1:n.islands]

    }

  }

  ### initialize the generations at which chromosomes will migrate between islands ###
  if (island.cluster.size != 1){

    migration.gens <- seq(migration.generations, generations - migration.generations, migration.generations)

  } else {

    migration.gens <- generations

  }


  ### evolve populations over island clusters ###
  first.seeds <- seq(1, n.islands, island.cluster.size)
  try(bplapply(first.seeds, function(cluster.start){
   # (lapply(first.seeds, function(cluster.start){

   cluster.seeds <- starting.seeds[cluster.start:(cluster.start + island.cluster.size - 1)]

   if (island.cluster.size > 1){

     ## initialize populations in the island cluster ##
     island.populations <- lapply(1:length(cluster.seeds), function(cluster.idx){

       cluster.seed.val <- cluster.seeds[cluster.idx]
       evolve.island(n.migrations = n.migrations, case.genetic.data = case.genetic.data,
                     complement.genetic.data = complement.genetic.data,
                     case.comp.different = case.comp.different,
                     case.minus.comp = case.minus.comp, both.one.mat = both.one.mat,
                     chrom.mat = chrom.mat, n.chromosomes = n.chromosomes,
                     n.candidate.snps = ncol(case.genetic.data), chromosome.size = chromosome.size,
                     start.generation = 1,
                     seed.val = cluster.seed.val, snp.chisq = snp.chisq,
                     original.col.numbers = original.col.numbers,
                     n.different.snps.weight = n.different.snps.weight,
                     n.both.one.weight = n.both.one.weight,
                     weight.function = weight.function, total.generations = migration.generations,
                     gen.same.fitness = gen.same.fitness,
                     max.generations = generations,
                     tol = tol, n.top.chroms = n.top.chroms, initial.sample.duplicates = initial.sample.duplicates,
                     snp.sampling.type = snp.sampling.type, crossover.prop = crossover.prop)

     })

     ### if we do not get convergence for all islands, migrate chromosomes ###
     all.converged <- F
     max.generations <- F
     while(!max.generations){

       for (island in 1:island.cluster.size){

         if (island == 1){

           island.populations[[island]]$chromosome.list <- c(island.populations[[island]]$chromosome.list,
                                                             island.populations[[island.cluster.size]]$migrations)

         } else {

           island.populations[[island]]$chromosome.list <- c(island.populations[[island]]$chromosome.list,
                                                             island.populations[[island - 1]]$migrations)

         }
       }

       ## evolving islands using the existing populations ##
       island.populations <- lapply(1:length(cluster.seeds), function(cluster.idx){

         cluster.seed.val <- cluster.seeds[cluster.idx]
         island <- island.populations[[cluster.idx]]
         evolve.island(n.migrations = n.migrations, case.genetic.data = case.genetic.data,
                       complement.genetic.data = complement.genetic.data,
                       case.comp.different = case.comp.different,
                       case.minus.comp = case.minus.comp, both.one.mat = both.one.mat,
                       chrom.mat = chrom.mat, n.chromosomes = n.chromosomes,
                       n.candidate.snps = ncol(case.genetic.data), chromosome.size = chromosome.size,
                       start.generation = island$generation,
                       seed.val = cluster.seed.val, snp.chisq = snp.chisq,
                       original.col.numbers = original.col.numbers, all.converged = all.converged,
                       n.different.snps.weight = n.different.snps.weight,
                       n.both.one.weight = n.both.one.weight,
                       weight.function = weight.function, total.generations = migration.generations,
                       gen.same.fitness = gen.same.fitness,
                       max.generations = generations,
                       tol = tol, n.top.chroms = n.top.chroms, initial.sample.duplicates = initial.sample.duplicates,
                       snp.sampling.type = snp.sampling.type, crossover.prop = crossover.prop,
                       chromosome.list = island$chromosome.list, fitness.score.mat = island$fitness.score.mat,
                       top.fitness = island$top.fitness, last.gens.equal = island$last.gens.equal,
                       top.generation.chromosome = island$top.generation.chromosome,
                       chromosome.mat.list = island$chromosome.mat.list,
                       sum.dif.vec.list = island$sum.dif.vec.list)

       })
       all.converged <- all(unlist(lapply(island.populations, function(x) x$last.gens.equal)))
       max.generations <- "top.chromosome.results" %in% names(island.populations[[1]])

     }

   } else {

     island.populations <- evolve.island(n.migrations = NULL, case.genetic.data = case.genetic.data,
                               complement.genetic.data = complement.genetic.data,
                               case.comp.different = case.comp.different,
                               case.minus.comp = case.minus.comp, both.one.mat = both.one.mat,
                               chrom.mat = chrom.mat, n.chromosomes = n.chromosomes,
                               n.candidate.snps = ncol(case.genetic.data), chromosome.size = chromosome.size,
                               start.generation = 1,
                               seed.val = cluster.seeds, snp.chisq = snp.chisq,
                               original.col.numbers = original.col.numbers,
                               n.different.snps.weight = n.different.snps.weight,
                               n.both.one.weight = n.both.one.weight,
                               weight.function = weight.function, total.generations = generations,
                               gen.same.fitness = gen.same.fitness,
                               max.generations = generations,
                               tol = tol, n.top.chroms = n.top.chroms, initial.sample.duplicates = initial.sample.duplicates,
                               snp.sampling.type = snp.sampling.type, crossover.prop = crossover.prop)

   }


    ### write results to file
    lapply(1:length(cluster.seeds), function(cluster.idx){

      cluster.seed.val <- cluster.seeds[cluster.idx]
      island <- island.populations[[cluster.idx]]
      out.file <- file.path(results.dir, paste0("island", cluster.seed.val, ".rds"))
      saveRDS(island, out.file)

    })

  }))
}


