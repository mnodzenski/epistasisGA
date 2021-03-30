#' A function to run the GADGETS method
#'
#' This function runs the GADGETS method on a given cluster of islands. It is a wrapper for
#' the underlying Rcpp function run_GADGETS.
#'
#' @param cluster.number An integer indicating the cluster number (used for labeling the output file).
#' @param results.dir The directory to which island results will be saved.
#' @param case.genetic.data The genetic data of the disease affected children from case-parent trios or affected/unaffected sibling pairs. Columns are SNP allele counts, and rows are individuals.
#' The ordering of the columns must be consistent with the LD structure specified in \code{block.ld.mat}.
#' @param complement.genetic.data A genetic dataset from the complements of the cases, where
#' \code{complement.genetic.data} = mother SNP counts + father SNP counts - case SNP counts.
#' Columns are SNP allele counts, rows are families. If using affected/unaffected sibling pairs, this should contain
#' the unaffected sibling genotypes.
#' @param case.comp.different A data frame or matrix indicating \code{case.genetic.data} != \code{complement.genetic.data},
#' where rows correspond to individuals and columns correspond to snps.
#' @param case.minus.comp A matrix equal to \code{case.genetic.data} - \code{complement genetic data}.
#' @param both.one.mat A matrix whose elements indicate whether both the case and complement have one copy of the minor allele,
#' equal to \code{case.genetic.data == 1 & complement.genetic.data == 1}.
#' @param block.ld.mat A logical, block diagonal matrix indicating whether the SNPs in \code{case.genetic.data} should be considered
#'  to be in linkage disequilibrium. Note that this means the ordering of the columns (SNPs) in \code{case.genetic.data} must be consistent
#'  with the LD blocks specified in \code{ld.block.mat}. In the absence of outside information, a reasonable default is to consider SNPs
#'  to be in LD if they are located on the same biological chromosome.
#' @param n.chromosomes An integer specifying the number of chromosomes to use in the GA.
#' @param chromosome.size An integer specifying the number of SNPs on each chromosome.
#' @param snp.chisq A vector of statistics to be used in sampling SNPs for mutation. By default, these are the square roots of
#' the chi-square marginal SNP-disease association statistics for each column in \code{case.genetic.data}, but can also be manually
#' specified or uniformly 1 (corresponding to totally random sampling).
#' @param original.col.numbers A vector of integers indicating the original column number of each SNP in \code{case.genetic.data}.
#' This is needed due to removal of low frequency SNPs in \code{preprocess.genetic.data}.
#' @param weight.lookup A vector that maps a family weight to the weighted sum of the number of different SNPs and SNPs both equal to one.
#' @param case2.mat A logical matrix indicating whether, for each SNP, the case carries 2 copies of the minor allele.
#' @param case0.mat A logical matrix indicating whether, for each SNP, the case carries 0 copies of the minor allele.
#' @param comp2.mat A logical matrix indicating whether, for each SNP, the complement/unaffected sibling carries 2 copies of the minor allele.
#' @param comp0.mat A logical matrix indicating whether, for each SNP, the complement/unaffected sibling carries 0 copies of the minor allele.
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
#' @param dif.coding A logical indicating whether, for a given SNP, the case - complement genotype difference should
#' be coded as the sign of the difference (defaulting to false) or the raw difference.
#' @return For each island in the cluster, an rds object containing a list with the following elements will be written to \code{results.dir}:
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
#' data(dad)
#' data(mom)
#' library(Matrix)
#' block.ld.mat <- as.matrix(bdiag(list(matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25))))
#' data.list <- preprocess.genetic.data(case[, 1:10], father.genetic.data = dad[ , 1:10],
#'                                mother.genetic.data = mom[ , 1:10],
#'                                block.ld.mat = block.ld.mat[1:10, 1:10])
#'
#'  case.genetic.data <- as.matrix(data.list$case.genetic.data)
#'  complement.genetic.data <- as.matrix(data.list$complement.genetic.data)
#'  original.col.numbers <- data.list$original.col.numbers
#'  chisq.stats <- data.list$chisq.stats
#'  block.ld.mat <- data.list$block.ld.mat
#'  case.minus.comp <- sign(as.matrix(case.genetic.data - complement.genetic.data))
#'  case.comp.different <- case.minus.comp != 0
#'  both.one.mat <- complement.genetic.data == 1 & case.genetic.data == 1
#'  case2.mat <- case.genetic.data == 2
#'  case0.mat <- case.genetic.data == 0
#'  comp2.mat <- complement.genetic.data == 2
#'  comp0.mat <- complement.genetic.data == 0
#'  snp.chisq <- sqrt(chisq.stats)
#'  weight.lookup <- vapply(seq_len(6), function(x) 2^x, 1)
#'  dir.create('tmp')
#'  GADGETS(cluster.number = 1, results.dir = 'tmp', case.genetic.data = case.genetic.data,
#'                    complement.genetic.data = complement.genetic.data,
#'                    case.comp.different = case.comp.different,
#'                    case.minus.comp = case.minus.comp, both.one.mat = both.one.mat,
#'                    block.ld.mat = block.ld.mat, n.chromosomes = 10,
#'                    chromosome.size = 3, snp.chisq = snp.chisq,
#'                    original.col.numbers = original.col.numbers,
#'                    weight.lookup = weight.lookup, case2.mat = case2.mat,
#'                    case0.mat = case0.mat, comp2.mat = comp2.mat,
#'                    comp0.mat = comp0.mat, n.migrations = 2,
#'                    migration.interval = 5, max.generations = 10)
#' unlink('tmp', recursive = TRUE)
#'
#' @importFrom data.table as.data.table setorder setDT rbindlist transpose
#' @useDynLib epistasisGA
#' @export

GADGETS <- function(cluster.number, results.dir , case.genetic.data, complement.genetic.data, case.comp.different,
                   case.minus.comp, both.one.mat, block.ld.mat, n.chromosomes, chromosome.size,
                   snp.chisq, original.col.numbers, weight.lookup, case2.mat, case0.mat, comp2.mat,
                   comp0.mat, island.cluster.size = 4, n.migrations = 20, n.different.snps.weight = 2,
                   n.both.one.weight = 1, migration.interval = 50, gen.same.fitness = 50,
                   max.generations = 500, initial.sample.duplicates = FALSE,
                   crossover.prop = 0.8, recessive.ref.prop = 0.75, recode.test.stat = 1.64,
                   dif.coding = FALSE) {

    ### run rcpp version of GADGET ##
    rcpp.res <- run_GADGETS(island.cluster.size, n.migrations, case.genetic.data,
                           complement.genetic.data, case.comp.different, case.minus.comp,
                           both.one.mat, block.ld.mat, n.chromosomes, chromosome.size,
                           weight.lookup, case2.mat, case0.mat, comp2.mat, comp0.mat,
                           snp.chisq, original.col.numbers,
                           n.different.snps.weight, n.both.one.weight, migration.interval,
                           gen.same.fitness, max.generations,
                           initial.sample.duplicates, crossover.prop, recessive.ref.prop,
                           recode.test.stat, dif.coding)

    ### clean up and output results
    lapply(seq_along(rcpp.res), function(island.number){

        #pick out the pieces from rcpp output
        n.generations <- rcpp.res[[island.number]][["generation"]]
        final.population.list <- rcpp.res[[island.number]][["current_fitness"]]
        chromosome.list <- final.population.list[["gen_original_cols"]]
        chromosome.dt <- as.data.table(do.call(rbind, chromosome.list))
        colnames(chromosome.dt) <- paste0("snp", seq_len(chromosome.size))
        fitness.score.dt <- data.table(fitness.score = final.population.list[["fitness_scores"]])
        dif.vec.list <- final.population.list[["sum_dif_vecs"]]
        dif.vec.dt <- as.data.table(do.call(rbind, dif.vec.list))
        colnames(dif.vec.dt) <- paste0("snp", seq_len(chromosome.size), ".diff.vec")
        risk.allele.vec.list <- final.population.list[["risk_allele_vecs"]]
        risk.allele.vec.dt <- as.data.table(do.call(rbind, risk.allele.vec.list))
        colnames(risk.allele.vec.dt) <- paste0("snp", seq_len(chromosome.size), ".allele.copies")
        n.case.risk.geno.dt <- data.table(n.cases.risk.geno = final.population.list[["n_case_risk_geno_vec"]])
        n.comp.risk.geno.dt <- data.table(n.comps.risk.geno = final.population.list[["n_comp_risk_geno_vec"]])
        final.result <- cbind(chromosome.dt, dif.vec.dt, risk.allele.vec.dt, fitness.score.dt,
                              n.case.risk.geno.dt, n.comp.risk.geno.dt)
        setorder(final.result, -fitness.score)

        #output list
        final.list <- list(top.chromosome.results = final.result, n.generations = n.generations)

        #write to file
        out.file <- file.path(results.dir, paste0("cluster", cluster.number, ".island", island.number,".rds"))
        saveRDS(final.list, out.file)

    })

}


