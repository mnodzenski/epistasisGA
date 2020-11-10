#' A function to run the GADGET method
#'
#' This function runs the GADGET method on a given cluster of islands. It is a wrapper for
#' the underlying C++ function run_GADGET.
#'
#' @param cluster.number An integer indicating the cluster number (used for labeling the output file).
#' @param results.dir The directory to which island results will be saved.
#' @param case.genetic.data The genetic data of the disease affected children from case-parent trios. Columns are SNPs, and rows are individuals.
#' The ordering of the columns must be consistent with the LD structure specified in \code{block.ld.mat}.
#' @param complement.genetic.data A genetic dataset from the complements of the cases, where
#' \code{complement.genetic.data} = mother SNP counts + father SNP counts - case SNP counts.
#' Columns are SNPs, rows are families.
#' @param case.comp.different A data frame or matrix indicating \code{case.genetic.data} != \code{complement.genetic.data},
#' where rows correspond to individuals and columns correspond to snps.
#' @param case.minus.comp A matrix equal to \code{case.genetic.data} - \code{complement genetic data}.
#' @param both.one.mat A matrix whose elements indicate whether both the case and complement have one copy of the minor allele,
#' equal to \code{case.genetic.data == 1 & complement.genetic.data == 1}.
#' @param block.ld.mat A logical, block diagonal matrix indicating whether the SNPs in \code{case.genetic.data} should be considered
#'  to be in linkage disequilibrium. Note that this means the ordering of the columns (SNPs) in \code{case.genetic.data} must be consistent
#'  with the LD blocks specified in \code{ld.block.mat}. In the absence of outside information, a reasonable default is to consider SNPs
#'  to be in LD if they are located on the same biological chromosome. If investigating maternal effects, where SNPs are being used as a
#'  proxy for a prenatal exposure, every entry of \code{block.ld.mat} should be set to TRUE.
#' @param n.chromosomes An integer specifying the number of chromosomes to use in the GA.
#' @param chromosome.size An integer specifying the number of SNPs on each chromosome.
#' @param snp.chisq A vector of chi-square statistics corresponding to marginal SNP-disease associations for each column in \code{case.genetic.data}.
#' @param original.col.numbers A vector of integers indicating the original column number of each SNP in \code{case.genetic.data},
#' needed due to removal of low frequency SNPs in \code{preprocess.genetic.data}.
#' @param weight.lookup A vector that maps a family weight to the weighted sum of the number of different SNPs and SNPs both equal to one.
#' @param island.cluster.size An integer specifying the number of islands in the cluster. See code{run.ga} for additional details.
#' @param n.migrations The number of chromosomes that migrate among islands. This value must be less than \code{n.chromosomes} and greater than 0, defaulting to 20.
#' @param n.different.snps.weight The number by which the number of different SNPs between a case and complement is multiplied in computing the family weights. Defaults to 2.
#' @param n.both.one.weight The number by which the number of SNPs equal to 1 in both the case and complement is multiplied in computing the family weights. Defaults to 1.
#' @param migration.interval The interval of generations for which the GA will run prior to migration of top chromosomes among islands in a cluster. Defaults to 50.
#' In other words, top chromosomes will migrate among cluster islands every \code{migration.interval} generations.
#' @param gen.same.fitness The number of consecutive generations with the same fitness score required for algorithm termination. Defaults to 50.
#' @param max.generations The maximum number of generations for which the GA will run. Defaults to 500.
#' @param tol The maximum absolute pairwise difference among the top fitness scores from the previous \code{gen.same.fitness} generations
#' considered to be sufficient to stop the algorithm.
#' @param n.top.chroms The number of top scoring chromosomes according to fitness score to return. Defaults to 100.
#' @param initial.sample.duplicates A logical indicating whether the same SNP can appear in more than one chromosome in the initial sample of chromosomes
#'  (the same SNP may appear in more than one chromosome thereafter, regardless). Default to FALSE.
#' @param crossover.prop A numeric between 0 and 1 indicating the proportion of chromosomes to be subjected to cross over.
#' The remaining proportion will be mutated. Defaults to 0.8.
#' @param n.case.high.risk.thresh The number of cases with the provisional high risk set required to check for recessive patterns of allele inheritance.
#' @param outlier.sd The number of standard deviations from the mean allele count used to determine whether recessive allele coding is appropriate
#' for a given SNP. See the GADGET paper for specific details on the implementation of this argument.
#' @return For each island in the cluster, an rds object containing a list with the following elements will be written to \code{results.dir}:
#' \describe{
#'  \item{top.chromosome.results}{A data.table of the top \code{n.top.chroms scoring chromosomes}, their fitness scores, their difference vectors,
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
#'  case.genetic.data <- data.list$case.genetic.data
#'  complement.genetic.data <- data.list$complement.genetic.data
#'  original.col.numbers <- data.list$original.col.numbers
#'  chisq.stats <- data.list$chisq.stats
#'  block.ld.mat <- data.list$block.ld.mat
#'  case.minus.comp <- sign(as.matrix(case.genetic.data - complement.genetic.data))
#'  case.comp.different <- case.minus.comp != 0
#'  both.one.mat <- complement.genetic.data == 1 & case.genetic.data == 1
#'  snp.chisq <- sqrt(chisq.stats)
#'  weight.lookup <- vapply(seq_len(6), function(x) 2^x, 1)
#'  dir.create('tmp')
#'  GADGET(cluster.number = 1, results.dir = 'tmp', case.genetic.data = case.genetic.data,
#'                    complement.genetic.data = complement.genetic.data,
#'                    case.comp.different = case.comp.different,
#'                    case.minus.comp = case.minus.comp, both.one.mat = both.one.mat,
#'                    block.ld.mat = block.ld.mat, n.chromosomes = 10,
#'                    chromosome.size = 3, snp.chisq = snp.chisq,
#'                    original.col.numbers = original.col.numbers,
#'                    weight.lookup = weight.lookup, n.migrations = 2,
#'                    migration.interval = 5, max.generations = 10)
#' unlink('tmp', recursive = TRUE)
#'
#' @importFrom data.table as.data.table setorder setDT rbindlist transpose
#' @useDynLib snpGADGET
#' @export

GADGET <- function(cluster.number, results.dir , case.genetic.data, complement.genetic.data, case.comp.different,
                   case.minus.comp, both.one.mat, block.ld.mat, n.chromosomes, chromosome.size,
                   snp.chisq, original.col.numbers, weight.lookup, island.cluster.size = 4, n.migrations = 20,
                   n.different.snps.weight = 2, n.both.one.weight = 1, migration.interval = 50, gen.same.fitness = 50,
                   max.generations = 500, tol = 10^-6, n.top.chroms = 100, initial.sample.duplicates = FALSE,
                   crossover.prop = 0.8, n.case.high.risk.thresh = 20, outlier.sd = 2.5) {

    ### run rcpp version of GADGET ##
    rcpp.res <- run_GADGET(island.cluster.size, n.migrations, case.genetic.data,
                           complement.genetic.data, case.comp.different, case.minus.comp,
                           both.one.mat, block.ld.mat, n.chromosomes, chromosome.size,
                           weight.lookup, snp.chisq, original.col.numbers,
                           n.different.snps.weight, n.both.one.weight, migration.interval,
                           gen.same.fitness, max.generations, tol, n.top.chroms,
                           initial.sample.duplicates, crossover.prop, n.case.high.risk.thresh,
                           outlier.sd)

    ### clean up and output results
    lapply(seq_along(rcpp.res), function(island.number){

        #pick out the pieces from rcpp output
        n.generations <- rcpp.res[[island.number]][["generation"]]
        fitness.score.vec <- unlist(rcpp.res[[island.number]][["fitness_score_list"]][seq_len(n.generations)])
        all.chrom.dt <- rbindlist(lapply(rcpp.res[[island.number]][["gen_chromosome_list"]][seq_len(n.generations)],
                                    function(gen.list) transpose(setDT(gen.list))))
        sum.dif.vec.dt <- rbindlist(lapply(rcpp.res[[island.number]][["sum_dif_vec_list"]][seq_len(n.generations)],
                                        function(gen.list) transpose(setDT(gen.list))))
        risk.allele.dt <- rbindlist(lapply(rcpp.res[[island.number]][["risk_allele_vec_list"]][seq_len(n.generations)],
                                        function(gen.list) transpose(setDT(gen.list))))

        unique.chromosome.dt <- unique(all.chrom.dt)
        colnames(unique.chromosome.dt) <- paste0("snp", seq_len(ncol(unique.chromosome.dt)))
        unique.chrom.dif.vec.dt <- sum.dif.vec.dt[!duplicated(all.chrom.dt), ]
        unique.chrom.risk.allele.vec.dt <-  risk.allele.dt[!duplicated(all.chrom.dt), ]
        colnames(unique.chrom.dif.vec.dt) <- paste0("snp", seq_len(ncol(unique.chrom.dif.vec.dt)),
                                                    ".diff.vec")
        colnames(unique.chrom.risk.allele.vec.dt) <- paste0("snp", seq_len(ncol(unique.chrom.dif.vec.dt)),
                                                            ".allele.copies")
        unique.fitness.score.vec <- fitness.score.vec[!duplicated(all.chrom.dt)]
        unique.results <- cbind(unique.chromosome.dt, unique.chrom.dif.vec.dt, unique.chrom.risk.allele.vec.dt)
        unique.results[, `:=`(raw.fitness.score, unique.fitness.score.vec)]
        unique.results[, `:=`(min.elem, min(abs(.SD))), by = seq_len(nrow(unique.results)),
                       .SDcols = (1 + chromosome.size):(2 * chromosome.size)]
        unique.results[, `:=`(fitness.score, min.elem * raw.fitness.score)]
        setorder(unique.results, -fitness.score)
        final.result <- unique.results[seq_len(n.top.chroms), ]

        #output list
        final.list <- list(top.chromosome.results = final.result, n.generations = n.generations)

        #write to file
        out.file <- file.path(results.dir, paste0("cluster", cluster.number, ".island", island.number,".rds"))
        saveRDS(final.list, out.file)

    })

}


