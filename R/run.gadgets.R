#' A function to run the GADGETS algorithm to detect multi-SNP effects in case-parent triad studies.
#'
#' This function runs the GADGETS algorithm to detect multi-SNP effects in case-parent triad studies.
#'
#' @param data.list The output list from \code{preprocess.genetic.data}.
#' @param n.chromosomes An integer specifying the number of chromosomes to use for each island in GADGETS.
#' @param chromosome.size An integer specifying the number of SNPs in each chromosome.
#' @param results.dir The directory to which island results will be saved.
#' @param cluster.type A character string indicating the type of cluster on which to evolve solutions in parallel.
#' Supported options are interactive, socket, multicore, sge, slurm, lsf, openlava, or torque. See the documentation for package batchtools for more information.
#' @param registryargs A list of the arguments to be provided to \code{batchtools::makeRegistry}.
#' @param resources A named list of key-value pairs to be substituted into the template file. Options available are specified in \code{batchtools::submitJobs}.
#' @param cluster.template A character string of the path to the template file required for the cluster specified in \code{cluster.type}.
#'  Defaults to NULL. Required for options sge, slurm, lsf, openlava and torque of argument \code{cluster.type}.
#' @param n.workers An integer indicating the number of workers for the cluster specified in \code{cluster.type}, if socket or multicore.
#' Defaults to \code{parallel::detectCores - 2}.
#' @param n.chunks An integer specifying the number of chunks jobs running island clusters should be split into when dispatching jobs using \code{batchtools}.
#' For multicore or socket \code{cluster.type}, this defaults to \code{n.workers}, resulting in the total number of island cluster jobs
#' (equal to \code{n.islands}\\\code{island.cluster.size}) being split into \code{n.chunks} chunks.
#' All chunks then run in parallel, with jobs within a chunk running sequentially. For other cluster types, this defaults to 1 chunk, with the recommendation
#' that users of HPC clusters which support array jobs specify \code{chunks.as.arrayjobs = TRUE} in argument \code{resources}. For those users, the setup will
#' submit an array of \code{n.islands}\\\code{island.cluster.size} jobs to the cluster. For HPC clusters that do not support array jobs, the default setting
#' should not be used. See \code{batchtools::submitJobs} for more information on job chunking.
#' @param n.different.snps.weight The number by which the number of different SNPs between a case and complement or unaffected sibling
#'  is multiplied in computing the family weights. Defaults to 2.
#' @param n.both.one.weight The number by which the number of SNPs equal to 1 in both the case and complement or unaffected sibling
#'  is multiplied in computing the family weights. Defaults to 1.
#' @param weight.function.int An integer used to assign family weights. Specifically, we use \code{weight.function.int} in a  function that takes the weighted sum
#' of the number of different SNPs and SNPs both equal to one as an argument, denoted as x, and returns a family weight equal to \code{weight.function.int}^x. Defaults to 2.
#' If set to null, then the family weight will not be exponentiated and instead set to just x.
#' @param generations The maximum number of generations for which GADGETS will run. Defaults to 500.
#' @param gen.same.fitness The number of consecutive generations with the same fitness score required for algorithm termination. Defaults to 50.
#' @param initial.sample.duplicates A logical indicating whether the same SNP can appear in more than one chromosome in the initial sample of chromosomes
#'  (the same SNP may appear in more than one chromosome thereafter, regardless). Default to FALSE.
#' @param snp.sampling.type A string indicating how SNPs are to be sampled for mutations. Options are 'chisq', 'random', or 'manual'. The 'chisq' option takes
#' into account the marginal association between a SNP and disease status, with larger marginal associations corresponding to higher sampling probabilities.
#' The 'random'  option gives each SNP the same sampling probability regardless of marginal association. The 'manual' option should be used when
#' \code{snp.sampling.probs} are manually input into function \code{preprocess.genetic.data}. Defaults to 'chisq'.
#' @param crossover.prop A numeric between 0 and 1 indicating the proportion of chromosomes to be subjected to cross over.
#' The remaining proportion will be mutated. Defaults to 0.8.
#' @param n.islands An integer indicating the number of islands to be used. Defaults to 1000.
#' @param island.cluster.size An integer specifying the number of islands in a given cluster. Must evenly divide \code{n.islands} and defaults to 4.
#' More specifically, under the default settings, the 1000 \code{n.islands} are split into 250 distinct clusters each containing 4 islands (\code{island.cluster.size}).
#' Within a cluster, migrations of top chromosomes from one cluster island to another are periodically permitted (controlled by \code{migration.generations}), and distinct
#' clusters evolve completely independently.
#' @param migration.generations An integer equal to the number of generations between migrations among islands of a distinct cluster.
#' Argument \code{generations} must be an integer multiple of this value. Defaults to 50.
#' @param n.migrations The number of chromosomes that migrate among islands. This value must be less than \code{n.chromosomes} and greater than 0, defaulting to 20.
#' @param recessive.ref.prop The proportion to which the observed proportion of informative cases with the provisional risk genotype(s) will be compared
#' to determine whether to recode the SNP as recessive. Defaults to 0.75.
#' @param recode.test.stat For a given SNP, the minimum test statistic required to recode and recompute the fitness score using recessive coding. Defaults to 1.64.
#' See the GADGETS paper for specific details.
#' @param use.parents A logical indicating whether parent data should be used in computing the fitness score. Defaults to FALSE. This should only be set to true
#' if the population is homogenous with no exposure related population structure.
#' @param use.parents.only A logical indicating whether only parent data should be used in computing the fitness score. Defaults to TRUE. This should only be set to true
#' if the population is homogenous with no exposure related population structure.
#' @return For each island, a list of two elements will be written to \code{results.dir}:
#' \describe{
#'  \item{top.chromosome.results}{A data.table of the final generation chromosomes, their fitness scores, their difference vectors,
#' and the number of risk alleles required for each chromosome SNP for a case or complement to be classified as having the provisional risk set.
#' See the package vignette for an example and the documentation for \code{chrom.fitness.score} for additional details.}
#'  \item{n.generations}{The total number of generations run.}
#' }
#'
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
#'                                ld.block.vec = c(10))
#' run.gadgets(pp.list, n.chromosomes = 4, chromosome.size = 3, results.dir = 'tmp',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'
#' unlink('tmp_bm', recursive = TRUE)
#' unlink('tmp', recursive = TRUE)
#' unlink('tmp_reg', recursive = TRUE)
#'
#' @importFrom bigmemory attach.big.matrix
#' @importFrom batchtools chunk makeRegistry batchMap submitJobs loadRegistry clearRegistry
#' @importFrom parallel detectCores
#' @export

run.gadgets <- function(data.list, n.chromosomes, chromosome.size, results.dir, cluster.type, registryargs = list(file.dir = NA,
    seed = 1500), resources = list(), cluster.template = NULL, n.workers = min(detectCores() - 2, n.islands/island.cluster.size),
    n.chunks = NULL, n.different.snps.weight = 2, n.both.one.weight = 1, weight.function.int = 2,
    generations = 500, gen.same.fitness = 50, initial.sample.duplicates = FALSE,
    snp.sampling.type = "chisq", crossover.prop = 0.8, n.islands = 1000, island.cluster.size = 4, migration.generations = 50,
    n.migrations = 20, recessive.ref.prop = 0.75, recode.test.stat = 1.64, use.parents = FALSE, use.parents.only = TRUE) {

    ### make sure if island clusters exist, the migration interval is set properly ###
    if (island.cluster.size > 1 & migration.generations >= generations & island.cluster.size != 1) {

        stop("migration.generations must be less than generations. Specify island.cluster.size = 1 and n.migrations = 0 if no migrations are desired.")

    }
    if (n.migrations == 0 & island.cluster.size != 1) {

        stop("Specify island.cluster.size = 1 and n.migrations = 0 if no migrations are desired.")

    }

    if (n.migrations != 0 & island.cluster.size == 1) {

        stop("Specify island.cluster.size = 1 and n.migrations = 0 if no migrations are desired.")

    }

    if (migration.generations == 1) {

        stop("migration.generations must be greater than 1")

    }
    if (island.cluster.size > 1 & generations%%migration.generations != 0) {

        stop("generations must be an integer multiple of migration.generations.")

    }

    ### make sure the island cluster size divides the number of islands evenly ###
    if (island.cluster.size > 1 & n.islands%%island.cluster.size != 0) {

        stop("n.islands must be an integer multiple of island.cluster.size")

    }

    ### make sure number of migrations is specified properly ###
    if (island.cluster.size > 1 & n.migrations >= n.chromosomes) {

        stop("n.migrations must be less than n.chromosomes")

    }

    ### make sure the weight function integer is actually an integer ###
    if (as.integer(weight.function.int) != weight.function.int){

        stop("weight.function.int must be an integer")

    }

    ### if no migrations, correctly set the migration.interval
    if (n.migrations == 0){

        migration.generations <- generations
    }

    ### compute the weight lookup table ###
    max.sum <- max(n.different.snps.weight, n.both.one.weight)*chromosome.size
    if (!is.null(weight.function.int)){

        weight.lookup <- vapply(seq_len(max.sum), function(x) weight.function.int^x, 1)

    } else {

        weight.lookup <- seq_len(max.sum)

    }

    storage.mode(weight.lookup) <- "integer"

    ### set sampling type for mutation snps ###
    if (snp.sampling.type == "chisq") {

        snp.chisq <- sqrt(data.list$chisq.stats)

    } else if (snp.sampling.type == "random") {

        snp.chisq <- rep(1, length(data.list$chisq.stats))

    } else if (snp.sampling.type == "manual"){

        snp.chisq <- data.list$chisq.stats

    }

    ### determine if islands have already been evolved
    clusters <- seq(1, n.islands, by = island.cluster.size)
    n.clusters <- length(clusters)
    cluster.ids <- seq_len(n.clusters)
    clusters.to.run <- cluster.ids
    if (!dir.exists(results.dir)) {

        dir.create(results.dir, recursive = TRUE)

    } else {

        prev.islands <- list.files(results.dir, pattern = "cluster")
        prev.islands <- prev.islands[!prev.islands %in% c("Old", "old", "combined")]
        prev.clusters <- unique(as.numeric(gsub("cluster|.island[0-9].rds", "", prev.islands)))
        clusters.to.run <- setdiff(clusters.to.run, prev.clusters)
        n.clusters.to.run <- length(clusters.to.run)

        if (n.clusters.to.run == 0) {

            stop("All islands have already been evolved")

        }

    }

    ### evolve populations over island clusters ###

    # make registry for submitting batch jobs
    reg.dir <- file.path(registryargs$file.dir, "registry")
    reg.dir <- gsub("//", "/", reg.dir, fixed = TRUE)
    if (!dir.exists(reg.dir)){

        registry <- do.call(makeRegistry, registryargs)
        registry$cluster.functions <- switch(cluster.type, interactive = batchtools::makeClusterFunctionsInteractive(),
                                             socket = batchtools::makeClusterFunctionsSocket(n.workers),
                                             multicore = batchtools::makeClusterFunctionsMulticore(n.workers),
                                             sge = batchtools::makeClusterFunctionsSGE(template = cluster.template),
                                             slurm = batchtools::makeClusterFunctionsSlurm(template = cluster.template),
                                             lsf = batchtools::makeClusterFunctionsLSF(template = cluster.template),
                                             openlava = batchtools::makeClusterFunctionsOpenLava(template = cluster.template),
                                             torque = batchtools::makeClusterFunctionsTORQUE(template = cluster.template),
                                             default = stop("unsupported cluster type '", cluster, "'"))
    } else {

        warning(paste("Registry already exists, loading from", reg.dir))
        registry <- loadRegistry(reg.dir, writeable = TRUE)
        clearRegistry(reg = registry)

    }

    # specify number of chunks
    if (is.null(n.chunks)) {

        if (!cluster.type %in% c("socket", "multicore")) {

            n.chunks <- 1

        } else {

            n.chunks <- n.workers

        }

    }

    # write jobs to registry
    ids <- batchMap(GADGETS, cluster.number = cluster.ids, more.args = list(results.dir = results.dir, n.migrations = n.migrations,
        case.genetic.data = data.list$case.genetic.data, complement.genetic.data = data.list$complement.genetic.data,
        ld.block.vec = data.list$ld.block.vec, n.chromosomes = n.chromosomes, chromosome.size = chromosome.size, snp.chisq = snp.chisq,
        weight.lookup = weight.lookup, island.cluster.size = island.cluster.size, n.different.snps.weight = n.different.snps.weight,
        n.both.one.weight = n.both.one.weight, migration.interval = migration.generations, gen.same.fitness = gen.same.fitness,
        max.generations = generations, initial.sample.duplicates = initial.sample.duplicates, crossover.prop = crossover.prop,
        recessive.ref.prop = recessive.ref.prop, recode.test.stat = recode.test.stat, exposure.levels = data.list$exposure.levels,
        exposure.risk.levels = data.list$exposure.risk.levels, exposure = data.list$exposure, use.parents = use.parents,
        use.parents.only = use.parents.only),
        reg = registry)

    # chunk the jobs
    ids[, `:=`(chunk, chunk(job.id, n.chunks = n.chunks))]

    # submit the jobs, using variable cluster.to.run to restrict submitted jobs to those not previously
    # run
    submitJobs(ids = ids[clusters.to.run, ], reg = registry, resources = resources)
}


