#' A function to identify the best family weights for GADGETS and re-run islands that did not use those
#'
#' This function chooses the best family weights among those input for GADGETS and re-run islands that did not
#' use those best weights.
#'
#' @param results.dir The directory in which individual island results from \code{run.gadgets} are saved.
#' @param preprocessed.list The initial list produced by function \code{preprocess.genetic.data}.

#' @return A data.table containing the results aggregated across islands. Note these results be written to \code{results.dir}
#' as 'combined.island.unique.chromosome.results.rds'. See the package vignette for more detailed descriptions of the content
#' of each output column.
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' data(snp.annotations)
#'
#' pp.list <- preprocess.genetic.data(case[, 1:10], father.genetic.data = dad[ , 1:10],
#'                                mother.genetic.data = mom[ , 1:10],
#'                                ld.block.vec = c(10))
#'
#' run.gadgets(pp.list, n.chromosomes = 4, chromosome.size = 3, results.dir = 'tmp',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'
#' best.weights <- choose.best.weights('tmp', pp.list)
#'
#' unlink('tmp', recursive = TRUE)
#' unlink('tmp_reg', recursive = TRUE)
#'
#' @importFrom data.table rbindlist setkey setorder `:=` setDT
#' @export

choose.best.weights <- function(results.dir, preprocessed.list) {

    # list all islands in the results data
    island.names <- list.files(results.dir, pattern = "cluster", full.names = TRUE)

    # pick out top chrom from each island and compute pdt stat
    top.chrom.res <- lapply(island.names, function(island.file) {

        island <- gsub(".rds", "", basename(island.file))
        island.data <- readRDS(island.file)
        chrom.results <- island.data$top.chromosome.results

        #take top scorers
        chrom.results <- chrom.results[1, ]

        # compute pdt stat
        n.cases <- chrom.results$n.cases.risk.geno
        n.comps <- chrom.results$n.comps.risk.geno
        pdt <- (n.cases - n.comps)/sqrt(n.cases + n.comps)

        # get tuning parameter info
        family.weight.info <- island.data$family.weight.info
        additional.param.info <- island.data$additional.param.info
        sampled.params <- island.data$sampled.params
        return(list(pdt = pdt, family.weight.info = family.weight.info,
                    sampled.params = sampled.params,
                    additional.param.info = additional.param.info))

    })
    pdt.stats <- vapply(top.chrom.res, function(chrom.res){

        return(chrom.res$pdt)

    }, 1.0)
    weight.info <- lapply(top.chrom.res, function(chrom.res){

        return(chrom.res$family.weight.info)

    })
    top.pdt <- max(pdt.stats)
    top.pdt.chroms <- which(pdt.stats == top.pdt)
    chosen.top.chrom <- sample(top.pdt.chroms, 1)
    best.weights <- weight.info[[chosen.top.chrom]]
    best.weights.other.parms <- top.chrom.res[[chosen.top.chrom]]$additional.param.info
    best.weights.sampled.parms <- top.chrom.res[[chosen.top.chrom]]$sampled.params
    keep.these <- vapply(weight.info, function(x) all(x == best.weights), TRUE)

    # get rid of the islands that didn't use the best weights
    lapply(island.names[!keep.these], function(island.file){

        unlink(island.file)

    })

    # re-run deleted islands now using the better weights
    n.chromosomes <- best.weights.other.parms$n.chromosomes
    chromosome.size <- best.weights.other.parms$chromosome.size
    cluster.type <- best.weights.other.parms$cluster.type
    registryargs <- best.weights.other.parms$registryargs

    #change the seed for the new islands
    registryargs$seed <- registryargs$seed + 1

    resources <- best.weights.other.parms$resources
    cluster.template <- best.weights.other.parms$cluster.template
    n.workers <- best.weights.other.parms$n.workers
    n.chunks <- best.weights.other.parms$n.chunks
    n.different.snps.weight <- best.weights['n.different.snps.weight']
    n.both.one.weight <- best.weights['n.both.one.weight']
    weight.function.int <- best.weights['weight.function.int']
    generations <- best.weights.other.parms$max.generations
    migration.generations <- best.weights.other.parms$migration.interval
    gen.same.fitness <- best.weights.other.parms$gen.same.fitness
    initial.sample.duplicates <- best.weights.other.parms$initial.sample.duplicates
    snp.sampling.type <- best.weights.other.parms$snp.sampling.type
    crossover.prop <- best.weights.other.parms$crossover.prop
    n.islands <- best.weights.other.parms$n.islands
    island.cluster.size <- best.weights.other.parms$island.cluster.size
    n.migrations <- best.weights.other.parms$n.migrations
    recessive.ref.prop <- best.weights.other.parms$recessive.ref.prop
    recode.test.stat <- best.weights.other.parms$recode.test.stat

    run.gadgets(preprocessed.list, n.chromosomes, chromosome.size, results.dir, cluster.type, registryargs, resources,
                cluster.template, n.workers, n.chunks, n.different.snps.weight, n.both.one.weight, weight.function.int,
                generations, migration.generations, gen.same.fitness, initial.sample.duplicates, snp.sampling.type,
                crossover.prop, n.islands, island.cluster.size, n.migrations, recessive.ref.prop, recode.test.stat)

    # return the best weights
    return(list(family.weights = best.weights, sampled.params = best.weights.sampled.parms,
                param.info = best.weights.other.parms))

}

