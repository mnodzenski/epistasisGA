#' A function to run the GADGETS method
#'
#' This function runs the GADGETS method on a given cluster of islands. It is a wrapper for
#' the underlying Rcpp function run_GADGETS.
#'
#' @param cluster.number An integer indicating the cluster number (used for labeling the output file).
#' @param results.dir The directory to which island results will be saved.
#' @param case.genetic.data The genetic data of the disease affected children from case-parent trios or affected/unaffected sibling pairs.
#' Columns are SNP allele counts, and rows are individuals. This object should be of class 'matrix'. The ordering of the columns must be consistent
#' with the LD structure specified in \code{ld.block.vec}. The genotypes cannot be dosages imputed with uncertainty. If any data are missing for a particular
#' family for a particular SNP, that SNP's genotype should be coded as -9 for the entire family, (\code{case.genetic.data} and
#' \code{father.genetic.data}/\code{mother.genetic.data} or \code{case.genetic.data} and \code{complement.genetic.data}).
#' @param complement.genetic.data A genetic dataset from the complements of the cases, where
#' \code{complement.genetic.data} = mother SNP counts + father SNP counts - case SNP counts. If using affected/unaffected siblings
#' this argument should be the genotypes for the unaffected siblings. This object should be of class 'matrix'. Columns are SNP allele counts, rows are
#' families. If not specified, \code{father.genetic.data} and \code{mother.genetic.data} must be specified.
#' The genotypes cannot be dosages imputed with uncertainty. If any data are missing for a particular family for a particular SNP, that SNP's genotype
#' should be coded as -9 for the entire family (\code{case.genetic.data} and \code{complement.genetic.data}).
#' @param ld.block.vec An integer vector specifying the linkage blocks of the input SNPs. As an example, for 100 candidate SNPs, suppose
#' we specify \code{ld.block.vec <- c(25, 75, 100)}. This vector indicates that the input genetic data has 3 distinct linkage blocks, with
#' SNPs 1-25 in the first linkage block, 26-75 in the second block, and 76-100 in the third block. Note that this means the ordering of the columns (SNPs)
#' in \code{case.genetic.data} must be consistent with the LD blocks specified in \code{ld.block.vec}. In the absence of outside information,
#' a reasonable default is to consider SNPs to be in LD if they are located on the same biological chromosome.
#' @param n.chromosomes An integer specifying the number of chromosomes to use in the GA.
#' @param chromosome.size An integer specifying the number of SNPs on each chromosome.
#' @param snp.chisq A vector of statistics to be used in sampling SNPs for mutation. By default, these are the square roots of
#' the chi-square marginal SNP-disease association statistics for each column in \code{case.genetic.data}, but can also be manually
#' specified or uniformly 1 (corresponding to totally random sampling).
#' @param weight.lookup A vector that maps a family weight to the weighted sum of the number of different SNPs and SNPs both equal to one.
#' @param null.mean.vec A vector of estimated null means for each of the three components of the Mahalanobis distance based
#' GxE fitness score. For all other uses, this should be specified as rep(0, 3) and will not be used.
#' @param null.se.vec An estimated vector of null standard errors for the three components of the
#' GxE fitness score. For all other uses, this should be specified as rep(1, 3) and will not be used.
#' @param island.cluster.size An integer specifying the number of islands in the cluster. See code{run.gadgets} for additional details.
#' @param n.migrations The number of chromosomes that migrate among islands. This value must be less than \code{n.chromosomes} and greater than 0, defaulting to 20.
#' @param n.different.snps.weight The number by which the number of different SNPs between a case and complement is multiplied in computing the family weights. Defaults to 2.
#' @param n.both.one.weight The number by which the number of SNPs equal to 1 in both the case and complement is multiplied in computing the family weights. Defaults to 1.
#' @param migration.interval The interval of generations for which GADGETS will run prior to migration of top chromosomes among islands in a cluster. Defaults to 50.
#' In other words, top chromosomes will migrate among cluster islands every \code{migration.interval} generations. We also check for convergence at each of these intervals.
#' @param gen.same.fitness The number of consecutive generations with the same fitness score required for algorithm termination. Defaults to 50.
#' @param max.generations The maximum number of generations for which GADGETS will run. Defaults to 500.
#' @param initial.sample.duplicates A logical indicating whether the same SNP can appear in more than one chromosome in the initial sample of chromosomes
#'  (the same SNP may appear in more than one chromosome thereafter, regardless). Defaults to FALSE.
#' @param crossover.prop A numeric between 0 and 1 indicating the proportion of chromosomes to be subjected to cross over.
#' The remaining proportion will be mutated. Defaults to 0.8.
#' @param recessive.ref.prop The proportion to which the observed proportion of informative cases with the provisional risk genotype(s) will be compared
#' to determine whether to recode the SNP as recessive. Defaults to 0.75.
#' @param recode.test.stat For a given SNP, the minimum test statistic required to recode and recompute the fitness score using recessive coding. Defaults to 1.64.
#' See the GADGETS paper for specific details.
#' @param exposure.levels An integer vector corresponding to the unique environmental exposure categories
#' from \code{exposure}.
#' @param exposure An integer vector corresponding to environmental exposures of the cases.
#' @return For each island in the cluster, an rds object containing a list with the following elements will be written to \code{results.dir}.
#' @param use.parents A integer indicating whether family level informativeness should be used alongside transmissions in computing GxE fitness scores. Defaults to 1,
#' indicating family level informativeness will be used. Specify 0 to only use transmission data.
#' \describe{
#'  \item{top.chromosome.results}{A data.table of the final generation chromosomes, their fitness scores, their difference vectors,
#' and the number of risk alleles required for each chromosome SNP for a case or complement to be classified as having the provisional risk set.
#' See the package vignette for an example and the documentation for \code{chrom.fitness.score} for additional details.}
#'  \item{n.generations}{The total number of generations run.}
#' }
#'
#' @examples
#'
#' set.seed(10)
#' data(case)
#' case <- as.matrix(case)
#' data(dad)
#' dad <- as.matrix(dad)
#' data(mom)
#' mom <- as.matrix(mom)
#' data.list <- preprocess.genetic.data(case[, 1:10], father.genetic.data = dad[ , 1:10],
#'                                mother.genetic.data = mom[ , 1:10],
#'                                ld.block.vec = c(10))
#'
#'  chisq.stats <- sqrt(data.list$chisq.stats)
#'  ld.block.vec <- data.list$ld.block.vec
#'  case.genetic.data <- data.list$case.genetic.data
#'  complement.genetic.data <- data.list$complement.genetic.data
#'
#'  #required inputs but not actually used in function below
#'  case.genetic.data.n <- matrix(0.0, 1, 1)
#'  complement.genetic.data.n <- matrix(0.0, 1, 1)
#'  exposure.mat <- data.list$exposure.mat + 0.0
#'
#'  weight.lookup <- vapply(seq_len(6), function(x) 2^x, 1)
#'  dir.create('tmp')
#' GADGETS(cluster.number = 1, results.dir = 'tmp', case.genetic.data = case.genetic.data,
#'        complement.genetic.data = complement.genetic.data, case.genetic.data.n = case.genetic.data.n,
#'        complement.genetic.data.n = complement.genetic.data.n, exposure.mat = exposure.mat,
#'        weight.lookup.n = weight.lookup + 0.0, mother.snps = NULL, child.snps = NULL,
#'        ld.block.vec = ld.block.vec,
#'        n.chromosomes = 10, chromosome.size = 3, snp.chisq = chisq.stats,
#'        weight.lookup = weight.lookup, n.migrations = 2, migration.interval = 5,
#'        gen.same.fitness = 10, max.generations = 10, null.mean.vec = rep(0, 3),
#'        null.se.vec = rep(1, 3))
#'
#' @importFrom data.table as.data.table setorder setDT rbindlist transpose
#' @useDynLib epistasisGAGE
#' @export

GADGETS <- function(cluster.number, results.dir, case.genetic.data, complement.genetic.data,
                    case.genetic.data.n, complement.genetic.data.n, exposure.mat, weight.lookup.n, ld.block.vec, n.chromosomes,
                    chromosome.size, snp.chisq, weight.lookup, null.mean.vec, null.se.vec, mother.snps, child.snps,
                    island.cluster.size = 4, n.migrations = 20, n.different.snps.weight = 2, n.both.one.weight = 1, migration.interval = 50,
                    gen.same.fitness = 50, max.generations = 500, initial.sample.duplicates = FALSE, crossover.prop = 0.8, recessive.ref.prop = 0.75,
                    recode.test.stat = 1.64, exposure.levels = NULL, exposure = NULL, use.parents = 1, cont.GxE = FALSE) {

    ### run rcpp version of GADGETS ##
    # rcpp.res <- run_GADGETS(island.cluster.size, n.migrations, ld.block.vec, n.chromosomes, chromosome.size,
    #                         weight.lookup,  snp.chisq, case.genetic.data, complement.genetic.data, null.mean.vec,
    #                         null.se.vec, exposure.levels, exposure, n.different.snps.weight, n.both.one.weight,
    #                         migration.interval, gen.same.fitness, max.generations, initial.sample.duplicates,
    #                         crossover.prop, recessive.ref.prop, recode.test.stat, use.parents)
    rcpp.res <- run_GADGETS(island.cluster.size, n.migrations, ld.block.vec, n.chromosomes, chromosome.size,
                            weight.lookup,  snp.chisq, case.genetic.data, complement.genetic.data, case.genetic.data.n,
                            complement.genetic.data.n, exposure.mat, weight.lookup.n, null.mean.vec,
                            null.se.vec, mother.snps, child.snps, exposure.levels, exposure, n.different.snps.weight,
                            n.both.one.weight, migration.interval, gen.same.fitness, max.generations, initial.sample.duplicates,
                            crossover.prop, recessive.ref.prop, recode.test.stat, use.parents, cont.GxE)

    ### clean up and output results
    lapply(seq_along(rcpp.res), function(island.number){

        #pick out the pieces from rcpp output
        rcpp.res.length <- length(rcpp.res[[island.number]])
        n.generations <- rcpp.res[[island.number]][["generation"]]
        final.population.list <- rcpp.res[[island.number]][["current_fitness"]]
        chromosome.list <- final.population.list[["gen_original_cols"]]
        chromosome.dt <- as.data.table(do.call(rbind, chromosome.list))
        colnames(chromosome.dt) <- paste0("snp", seq_len(chromosome.size))
        fitness.score.dt <- data.table(fitness.score = final.population.list[["fitness_scores"]])
        dif.vec.list <- final.population.list[["sum_dif_vecs"]]
        dif.vec.dt <- as.data.table(do.call(rbind, dif.vec.list))
        colnames(dif.vec.dt) <- paste0("snp", seq_len(chromosome.size), ".diff.vec")

        if (!cont.GxE){

            risk.allele.vec.list <- final.population.list[["risk_allele_vecs"]]
            risk.allele.vec.dt <- as.data.table(do.call(rbind, risk.allele.vec.list))
            colnames(risk.allele.vec.dt) <- paste0("snp", seq_len(chromosome.size), ".allele.copies")
            n.case.risk.geno.dt <- data.table(n.cases.risk.geno = final.population.list[["n_case_risk_geno_vec"]])
            n.comp.risk.geno.dt <- data.table(n.comps.risk.geno = final.population.list[["n_comp_risk_geno_vec"]])
            final.result <- cbind(chromosome.dt, dif.vec.dt, risk.allele.vec.dt, fitness.score.dt,
                                  n.case.risk.geno.dt, n.comp.risk.geno.dt)
            setorder(final.result, -fitness.score)
        } else {

            final.result <- cbind(chromosome.dt, dif.vec.dt, fitness.score.dt)
            setorder(final.result, -fitness.score)

        }

        # } else {
        #
        #     # grab exposure level info
        #     exposure.info <- final.population.list[["exposure_level_info"]]
        #
        #     # put together results for each chrom
        #     exposure.info.dt <- rbindlist(lapply(exposure.info, function(chrom.res){
        #
        #         exposure.info.list <- chrom.res[["score_by_exposure"]]
        #         dt.list <- lapply(seq_along(exposure.info.list), function(exposure.number){
        #
        #             # pick out results
        #             exposure.level <- exposure.levels[exposure.number]
        #             exposure.res <- exposure.info.list[[exposure.number]]
        #
        #             # specify column names
        #             colnames.start <- paste0("exposure", exposure.level)
        #             diff.vec.colnames <- paste0("snp", seq_len(chromosome.size), ".diff.vec")
        #             allele.copy.colnames <- paste0("snp", seq_len(chromosome.size), ".allele.copies")
        #             colnames.end <- c(diff.vec.colnames, allele.copy.colnames,
        #                               "n.cases.risk.geno", "n.comps.risk.geno")
        #
        #             # pick out results
        #             if ("no_informative_families" %in% names(exposure.res)){
        #
        #                 res <- data.table(matrix(NA, 1, length(colnames.end)))
        #                 res[ , 1] <- risk.order[exposure.number]
        #
        #             } else {
        #
        #                 exposure.res.target <- exposure.res[c("sum_dif_vecs", "risk_set_alleles", "n_case_risk_geno",
        #                                                       "n_comp_risk_geno")]
        #                 res <- data.table(t(data.frame(unlist(exposure.res.target))))
        #
        #             }
        #             colnames(res) <- paste(colnames.start, colnames.end)
        #             return(res)
        #
        #         })
        #
        #         combined.dt <- setDT(unlist(dt.list, recursive = FALSE), check.names = TRUE)[]
        #         return(combined.dt)
        #
        #     }))
        #
        #
        #     final.result <- cbind(chromosome.dt, dif.vec.dt, fitness.score.dt, exposure.info.dt)
        #     setorder(final.result, -fitness.score)
        # }

        #output list
        final.list <- list(top.chromosome.results = final.result, n.generations = n.generations)

        #write to file
        out.file <- file.path(results.dir, paste0("cluster", cluster.number, ".island", island.number,".rds"))
        saveRDS(final.list, out.file)

    })

}


