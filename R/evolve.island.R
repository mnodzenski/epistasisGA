#' A function to run a genetic algorithm for a specific island for a given number of generations
#'
#' This function runs a genetic algorithm for a specific island for a given number of generations.
#' This function should not be used independently, it is only intended for use in \code{run.ga}.
#'
#' @param n.migrations The number of chromosomes that migrate among islands. This value must be less than \code{n.chromosomes} and greater than 0, defaulting to 20.
#' @param case.genetic.data The genetic data of the disease affected children from case-parent trios. Columns are SNPs, and rows are individuals.
#' @param complement.genetic.data A genetic dataset from the complements of the cases, where
#' \code{complement.genetic.data} = mother SNP counts + father SNP counts - case SNP counts.
#' Columns are SNPs, rows are families.
#' @param case.comp.different A data frame or matrix indicating \code{case.genetic.data} != \code{complement.genetic.data},
#' where rows correspond to individuals and columns correspond to snps.
#' @param case.minus.comp A matrix equal to \code{case.genetic.data} - \code{complement genetic data}.
#' @param both.one.mat A matrix whose elements indicate whether both the case and complement have one copy of the minor allele,
#' equal to \code{case.genetic.data == 1 & complement.genetic.data == 1}.
#' @param chrom.mat A logical matrix indicating whether the SNPs in \code{case.genetic.data} are located on the same biological chromosome.
#' @param n.chromosomes An integer specifying the number of chromosomes to use in the GA.
#' @param n.candidate.snps A scalar indicating the number eligible SNPs in the input data, after filtering out low MAF SNPs.
#' @param chromosome.size An integer specifying the number of SNPs on each chromosome.
#' @param start.generation The generation at which this function should begin. If 1, a random set of chromosomes will be initialized.
#' Otherwise the argument \code{chromosome.list} will be used.
#' @param snp.chisq A vector of chi-square statistics corresponding to marginal SNP-disease associations for each column in \code{case.genetic.data}.
#' @param original.col.numbers A vector of integers indicating the original column number of each SNP in \code{case.genetic.data},
#' needed due to removal of low frequency SNPs in \code{preprocess.genetic.data}.
#' @param all.converged A logical indicating whether each island in the cluster has converged to a solution.
#' @param n.different.snps.weight The number by which the number of different SNPs between a case and complement is multiplied in computing the family weights. Defaults to 2.
#' @param n.both.one.weight The number by which the number of SNPs equal to 1 in both the case and complement is multiplied in computing the family weights. Defaults to 1.
#' @param weight.function A function that takes the weighted sum of the number of different SNPs and SNPs both equal to one as an argument, denoted as x,
#'  and returns a family weight. Defaults to 2^x.
#' @param migration.interval The interval of generations for which the GA will run prior to migration of top chromosomes among islands in a cluster. Defaults to 50.
#' In other words, top chromosomes will migrate among cluster islands every \code{migration.interval} generations.
#' @param max.generations The maximum number of generations for which the GA will run. Defaults to 500.
#' @param gen.same.fitness The number of consecutive generations with the same fitness score required for algorithm termination. Defaults to 50.
#' @param tol The maximum absolute pairwise difference among the top fitness scores from the previous \code{gen.same.fitness} generations
#' considered to be sufficient to stop the algorithm.
#' @param n.top.chroms The number of top scoring chromosomes according to fitness score to return. Defaults to 100.
#' @param initial.sample.duplicates A logical indicating whether the same SNP can appear in more than one chromosome in the initial sample of chromosomes
#'  (the same SNP may appear in more than one chromosome thereafter, regardless). Default to FALSE.
#' @param snp.sampling.type A string indicating how SNPs are to be sampled for mutations. Options are 'chisq' or 'random'. The 'chisq' option takes
#' into account the marginal association between a SNP and disease status, with larger marginal associations corresponding to higher sampling probabilities.
#' The 'random'  option gives each SNP the same sampling probability regardless of marginal association. Defaults to 'chisq'.
#' @param crossover.prop A numeric between 0 and 1 indicating the proportion of chromosomes to be subjected to cross over.
#' The remaining proportion will be mutated. Defaults to 0.8.
#' @param chromosome.list A list of input chromosomes which the genetic algorithm will use as a starting population.
#' @param fitness.score.mat A matrix of fitness scores from previous generations in island evolution.
#' @param top.fitness A vector of top fitness scores from previous generations.
#' @param last.gens.equal A logical indicating whether the last \code{gen.same.fitness} generations of the algorithm all produced the same top fitness score.
#' @param top.generation.chromosome A list of top chromosomes from previous generations.
#' @param chromosome.mat.list A list of matrices containing all chromosomes from previous generations.
#' @param sum.dif.vec.list A list of matrices containing the sum of differences vectors for chromosomes in previous generations.
#' @param n.case.high.risk.thresh The number of cases with the provisional high risk set required to check for recessive patterns of allele inheritance.
#' @return A list containing the following:
#' \describe{
#'  \item{migrations}{A list of chromosomes that will migrate to another island in the cluster.}
#'  \item{chromosome.list}{A list of chromosomes that will remain the island after migrations.}
#'  \item{fitness.score.mat}{A matrix of fitness scores for each chromosome across generations.}
#'  \item{top.fitness}{A vector of the top fitness scores across generations.}
#'  \item{last.gens.equal}{A logical indicating whether the last \code{gen.same.fitness} generations of the algorithm all produced the same top fitness score.}
#'  \item{top.generation.chromosome}{A list of the top scoring chromosome for the previous generations.}
#'  \item{chromosome.mat.list}{A list of all of the chromosomes examined in previous generations.}
#'  \item{sum.dif.vec.list}{A list of matrices containing the sum of the difference vectors for all examined chromosomes.}
#'  \item{generation}{An integer corresponding to the current generation of the genetic algorithm.}
#' }
#'
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' library(Matrix)
#' chrom.mat <- as.matrix(bdiag(list(matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25))))
#' data.list <- preprocess.genetic.data(case[, 1:10], father.genetic.data = dad[ , 1:10],
#'                                mother.genetic.data = mom[ , 1:10],
#'                                chrom.mat = chrom.mat[ , 1:10])
#'
#'  case.genetic.data <- data.list$case.genetic.data
#'  complement.genetic.data <- data.list$complement.genetic.data
#'  original.col.numbers <- data.list$original.col.numbers
#'  chisq.stats <- data.list$chisq.stats
#'  chrom.mat <- data.list$chrom.mat
#'  case.minus.comp <- sign(as.matrix(case.genetic.data - complement.genetic.data))
#'  case.comp.different <- case.minus.comp != 0
#'  both.one.mat <- complement.genetic.data == 1 & case.genetic.data == 1
#'  snp.chisq <- sqrt(chisq.stats)
#'  ex.island <- evolve.island(n.migrations = 2, case.genetic.data = case.genetic.data,
#'                    complement.genetic.data = complement.genetic.data,
#'                    case.comp.different = case.comp.different,
#'                    case.minus.comp = case.minus.comp, both.one.mat = both.one.mat,
#'                    chrom.mat = chrom.mat, n.chromosomes = 10,
#'                    n.candidate.snps = ncol(case.genetic.data),
#'                    chromosome.size = 3, start.generation = 1,
#'                    snp.chisq = snp.chisq,
#'                    original.col.numbers = original.col.numbers,
#'                    migration.interval = 5, max.generations = 10)
#'
#'
#' @importFrom matrixStats colSds rowMaxs
#' @importFrom data.table data.table rbindlist setorder
#' @importFrom stats rbinom sd
#' @export

evolve.island <- function(n.migrations = 20, case.genetic.data, complement.genetic.data, case.comp.different,
    case.minus.comp, both.one.mat, chrom.mat, n.chromosomes, n.candidate.snps, chromosome.size, start.generation,
    snp.chisq, original.col.numbers, all.converged = FALSE, n.different.snps.weight = 2, n.both.one.weight = 1,
    weight.function = function(x) 2^x, migration.interval = 50, gen.same.fitness = 50, max.generations = 500,
    tol = 10^-6, n.top.chroms = 100, initial.sample.duplicates = FALSE, snp.sampling.type = "chisq",
    crossover.prop = 0.8, chromosome.list = NULL, fitness.score.mat = NULL, top.fitness = NULL, last.gens.equal = NULL,
    top.generation.chromosome = NULL, chromosome.mat.list = NULL, sum.dif.vec.list = NULL, n.case.high.risk.thresh = 20) {

    ### initialize groups of candidate solutions if generation 1 ###
    generation <- start.generation
    generations <- min(start.generation + migration.interval - 1, max.generations)

    ### initiate chromosome list if this is the first generation ###
    if (generation == 1) {

        if ((ncol(case.genetic.data) < n.chromosomes * chromosome.size) & !initial.sample.duplicates) {

            print("Not enough SNPs present to allow for no initial sample duplicate SNPs, now allowing initial sample duplicate snps.")
            initial.sample.duplicates <- TRUE

        }
        chromosome.list <- vector(mode = "list", length = n.chromosomes)
        all.snps.idx <- seq_len(n.candidate.snps)
        for (i in seq_len(n.chromosomes)) {

            snp.idx <- sort(sample(all.snps.idx, chromosome.size, replace = FALSE))
            if (!initial.sample.duplicates) {

                all.snps.idx <- all.snps.idx[!all.snps.idx %in% snp.idx]

            }
            chromosome.list[[i]] <- snp.idx

        }

        fitness.score.mat <- matrix(rep(NA, max.generations * n.chromosomes), nrow = max.generations)
        top.fitness <- rep(0, max.generations)
        last.gens.equal <- FALSE
        top.generation.chromosome <- vector(mode = "list", length = max.generations)
        chromosome.mat.list <- vector(mode = "list", length = max.generations)
        sum.dif.vec.list <- vector(mode = "list", length = max.generations)

    }

    ### iterate over generations ###
    while (generation <= generations & !all.converged) {

        # print(paste('generation', generation))

        ### 1. compute the fitness score for each set of candidate snps ### print('Step 1/9')

        fitness.score.list <- lapply(seq_len(length(chromosome.list)), function(x) {

            chrom.fitness.score(case.genetic.data, complement.genetic.data, case.comp.different, chromosome.list[[x]],
                case.minus.comp, both.one.mat, chrom.mat, n.different.snps.weight, n.both.one.weight,
                weight.function, n.case.high.risk.thresh)

        })

        fitness.scores <- sapply(fitness.score.list, function(x) x$fitness.score)
        sum.dif.vecs <- t(sapply(fitness.score.list, function(x) x$sum.dif.vecs))

        # store the fitness scores, elements (snps) of the chromosomes, sum of the difference vectors
        fitness.score.mat[generation, ] <- fitness.scores
        chromosome.mat.list[[generation]] <- data.table(t(sapply(chromosome.list, function(x) original.col.numbers[x])))
        sum.dif.vec.list[[generation]] <- data.table(sum.dif.vecs)

        ### 2. identify the top scoring candidate solution(s) and fitness score ### print('Step 2/9')
        max.fitness <- max(fitness.scores)
        top.chromosome.idx <- which(fitness.scores == max.fitness)
        top.chromosomes <- chromosome.list[top.chromosome.idx]

        # if we have the same chromosome multiple times, just take one of them and put the duplicates in
        # the pool to be resampled
        if (length(top.chromosomes) > 1) {

            top.chromosome <- top.chromosomes[1]
            duplicate.top.chromosomes <- top.chromosome.idx[-1]

        } else {

            top.chromosome <- top.chromosomes
            duplicate.top.chromosomes <- NULL

        }

        ### 3. identify the lower scoring chromosomes ### print('Step 3/9')
        lower.chromosomes <- c(which(fitness.scores != max.fitness), duplicate.top.chromosomes)

        ### 4. Sample with replacement from the existing chromosomes ### print('Step 4/9') allow the top
        ### scoring chromosome to be sampled, but only sample from the unique chromosomes available
        sample.these <- !duplicated(chromosome.list)
        sampled.lower.idx <- sample(which(sample.these), length(lower.chromosomes), replace = TRUE,
            prob = fitness.scores[sample.these])
        sampled.lower.chromosomes <- chromosome.list[sampled.lower.idx]
        sampled.lower.dif.vecs <- sum.dif.vecs[sampled.lower.idx, ]
        sampled.lower.fitness.scores <- fitness.scores[sampled.lower.idx]

        ### 5. Determine whether each lower chromosome will be subject to mutation or crossing over ###
        ### print('Step 5/9')

        # only allowing cross-overs between distinct chromosomes (i.e, if a chromosome was sampled twice,
        # it can't cross over with itself)
        unique.lower.idx <- unique(sampled.lower.idx)

        cross.overs <- rep(FALSE, length(unique.lower.idx))

        if (round(length(unique.lower.idx) * crossover.prop)%%2 == 0 | length(unique.lower.idx) ==
            1) {

            cross.overs[sample(seq_len(length(unique.lower.idx)), size = round(length(unique.lower.idx) *
                crossover.prop))] <- TRUE

        } else {

            cross.overs[sample(seq_len(length(unique.lower.idx)), size = (round(length(unique.lower.idx) *
                crossover.prop) + 1))] <- TRUE

        }

        # those not getting crossover will be mutated
        mutations <- !cross.overs

        ### 6. Execute crossing over for the relevant chromosomes ### print('Step 6/9')
        if (sum(cross.overs) > 1) {

            cross.over.positions <- match(unique.lower.idx[cross.overs], sampled.lower.idx)
            cross.over.starts <- seq(1, length(cross.over.positions), by = 2)
            for (i in cross.over.starts) {

                # grab pair of chromosomes to cross over
                chrom1 <- sampled.lower.chromosomes[[cross.over.positions[i]]]
                chrom1.dif.vecs <- sampled.lower.dif.vecs[cross.over.positions[i], ]
                chrom1.fitness.score <- sampled.lower.fitness.scores[cross.over.positions[i]]

                chrom2 <- sampled.lower.chromosomes[[cross.over.positions[i + 1]]]
                chrom2.dif.vecs <- sampled.lower.dif.vecs[cross.over.positions[i + 1], ]
                chrom2.fitness.score <- sampled.lower.fitness.scores[cross.over.positions[i + 1]]

                # check for overlapping snps
                c1.c2.matching.snp.positions <- match(chrom1, chrom2)
                c1.c2.matching.snp.positions <- c1.c2.matching.snp.positions[!is.na(c1.c2.matching.snp.positions)]
                c1.c2.not.matching.snp.positions <- setdiff(seq_len(chromosome.size), c1.c2.matching.snp.positions)

                c2.c1.matching.snp.positions <- match(chrom2, chrom1)
                c2.c1.matching.snp.positions <- c2.c1.matching.snp.positions[!is.na(c2.c1.matching.snp.positions)]
                c2.c1.not.matching.snp.positions <- setdiff(seq_len(chromosome.size), c2.c1.matching.snp.positions)

                # for the non-matching snps, order by the decreasing magnitude of the difference vector in the
                # chromsome with the higher fitness score and order by the increasing magntidue of the difference
                # vector in the chromosome with the lower fitness score **ultimately will be used to substitute the
                # higher magnitude elements in the lower scoring chromosome for the lower magnitude elements in the
                # higher scoring chromosome**
                if (chrom1.fitness.score >= chrom2.fitness.score) {

                  c1.c2.not.matching.snp.positions <- c1.c2.not.matching.snp.positions[order(abs(chrom2.dif.vecs[c1.c2.not.matching.snp.positions]))]
                  c2.c1.not.matching.snp.positions <- c2.c1.not.matching.snp.positions[order(abs(chrom1.dif.vecs[c2.c1.not.matching.snp.positions]),
                    decreasing = TRUE)]

                } else {

                  c1.c2.not.matching.snp.positions <- c1.c2.not.matching.snp.positions[order(abs(chrom2.dif.vecs[c1.c2.not.matching.snp.positions]),
                    decreasing = TRUE)]
                  c2.c1.not.matching.snp.positions <- c2.c1.not.matching.snp.positions[order(abs(chrom1.dif.vecs[c2.c1.not.matching.snp.positions]))]

                }

                # order the chromosomes, first by the overlapping snps and the non-overlapping
                chrom2 <- chrom2[c(c1.c2.matching.snp.positions, c1.c2.not.matching.snp.positions)]
                chrom1 <- chrom1[c(c2.c1.matching.snp.positions, c2.c1.not.matching.snp.positions)]

                # determine how many snps could be crossed over and also make sure we don't simply swap chromosomes
                possible.cut.points <- sort(which(chrom1 != chrom2), decreasing = TRUE)
                if (length(possible.cut.points) == chromosome.size) {

                  n.possible.crosses <- chromosome.size - 1

                } else {

                  n.possible.crosses <- length(possible.cut.points)

                }

                # determine how many snps will actually be crossed over
                n.crosses <- sample.int(n.possible.crosses, 1)

                # pick out their positions
                cross.points <- c(seq_len(chromosome.size))[(chromosome.size - n.crosses + 1):chromosome.size]

                # exchange the high magnitude elements from the lower scoring chromosome with the low magnitude
                # elements from the high scoring chromsosome
                chrom1.cross <- chrom1
                chrom1.cross[cross.points] <- chrom2[cross.points]
                chrom2.cross <- chrom2
                chrom2.cross[cross.points] <- chrom1[cross.points]

                # replace in the chromosome list
                sampled.lower.chromosomes[[cross.over.positions[i]]] <- chrom1.cross
                sampled.lower.chromosomes[[cross.over.positions[i + 1]]] <- chrom2.cross

            }

        } else {

            cross.over.positions <- 0

        }

        ### 7. Mutate the chromosomes that were not crossed over ### print('Step 7/9')
        if (cross.over.positions[1] == 0) {

            mutation.positions <- (seq_len(length(sampled.lower.chromosomes)))

        } else {

            mutation.positions <- (seq_len(length(sampled.lower.chromosomes)))[-cross.over.positions]

        }
        snps.for.mutation <- sample(seq_len(n.candidate.snps), n.candidate.snps, prob = snp.chisq,
            replace = TRUE)
        for (i in mutation.positions) {

            # grab the chromosome and its difference vector
            target.chrom <- sampled.lower.chromosomes[[i]]
            target.dif.vec <- as.vector(t(sampled.lower.dif.vecs[i, ]))

            # sort the chromosome elements from lowest absolute difference vector to highest
            target.chrom <- target.chrom[order(abs(target.dif.vec))]

            # remove the chromosome's snps from the pool of available snps and sample new snps for the
            # mutations
            possible.snps.for.mutation <- snps.for.mutation[!snps.for.mutation %in% target.chrom]

            # determine which snps to mutate
            total.mutations <- min(max(1, rbinom(1, chromosome.size, 0.5)), length(unique(possible.snps.for.mutation)))
            mutate.these <- seq_len(total.mutations)

            # execute mutations
            mutated.snps <- rep(NA, total.mutations)
            for (j in seq_len(total.mutations)) {

                sampled.snp <- sample(possible.snps.for.mutation, 1)
                mutated.snps[j] <- sampled.snp
                possible.snps.for.mutation <- possible.snps.for.mutation[possible.snps.for.mutation !=
                  sampled.snp]

            }

            # substitute in mutations
            target.chrom[mutate.these] <- mutated.snps
            sampled.lower.chromosomes[[i]] <- target.chrom

        }

        ### 8. Combine into new population (i.e., the final collection of chromosomes for the next
        ### generation) print('Step 8/9')
        chromosome.list <- lapply(c(top.chromosome, sampled.lower.chromosomes), sort)

        ### 9.Increment Iterators print('Step 9/9')
        top.fitness[generation] <- max.fitness
        top.generation.chromosome[[generation]] <- original.col.numbers[top.chromosome[[1]]]
        # print(paste0('Max fitness score:', max.fitness)) print('Top Chromosome(s):')
        # print(original.col.numbers[top.chromosome[[1]]])
        if (generation >= gen.same.fitness) {

            # check to see if enough of the last generations have had the same top chromosome to terminate
            last.gens <- top.fitness[(generation - (gen.same.fitness - 1)):generation]
            last.gens.equal <- abs(max(last.gens) - min(last.gens)) < tol
            if (is.null(n.migrations) & last.gens.equal) {

                all.converged <- TRUE

            }

        }

        generation <- generation + 1

    }
    ### If the algorithm hasn't hit the max number of generations or converged, return a partial list of
    ### the results ###
    if (generation < max.generations & !all.converged & !is.null(n.migrations)) {

        # pick out the top fitness scores, and bottom fitness scores
        fitness.score.list <- lapply(seq_len(length(chromosome.list)), function(x) {

            chrom.fitness.score(case.genetic.data, complement.genetic.data, case.comp.different, chromosome.list[[x]],
                case.minus.comp, both.one.mat, chrom.mat, n.different.snps.weight, n.both.one.weight,
                weight.function)

        })
        fitness.scores <- sapply(fitness.score.list, function(x) x$fitness.score)
        chromosome.list <- chromosome.list[order(fitness.scores, decreasing = TRUE)]

        ### identify the chromosomes that will migrate to other islands ###
        migrations <- chromosome.list[seq_len(n.migrations)]

        ### remove the lowest scoring chromosomes in preparation for migration ###
        chromosome.list <- chromosome.list[seq_len((length(chromosome.list) - n.migrations))]

        ### return list of results ###
        return(list(migrations = migrations, chromosome.list = chromosome.list, fitness.score.mat = fitness.score.mat,
            top.fitness = top.fitness, last.gens.equal = last.gens.equal, top.generation.chromosome = top.generation.chromosome,
            chromosome.mat.list = chromosome.mat.list, sum.dif.vec.list = sum.dif.vec.list, generation = generation))

    } else {

        ### Otherwise Return the best chromosomes, their fitness scores, difference vectors and the number of
        ### generations ###
        last.generation <- generation - 1
        all.chrom.dt <- rbindlist(chromosome.mat.list[seq_len(last.generation)], use.names = FALSE)
        all.chrom.dif.vec.dt <- rbindlist(sum.dif.vec.list[seq_len(last.generation)], use.names = FALSE)
        unique.chromosome.dt <- unique(all.chrom.dt)
        colnames(unique.chromosome.dt) <- paste0("snp", seq_len(ncol(unique.chromosome.dt)))
        unique.chrom.dif.vec.dt <- all.chrom.dif.vec.dt[!duplicated(all.chrom.dt), ]
        colnames(unique.chrom.dif.vec.dt) <- paste0("snp", seq_len(ncol(unique.chrom.dif.vec.dt)),
            ".diff.vec")
        unique.fitness.score.vec <- as.vector(t(fitness.score.mat[seq_len(last.generation), ]))[!duplicated(all.chrom.dt)]
        unique.results <- cbind(unique.chromosome.dt, unique.chrom.dif.vec.dt)
        unique.results[, `:=`(raw.fitness.score, unique.fitness.score.vec)]
        unique.results[, `:=`(min.elem, min(abs(.SD))), by = seq_len(nrow(unique.results)), .SDcols = (1 +
            chromosome.size):(2 * chromosome.size)]
        unique.results[, `:=`(fitness.score, min.elem * raw.fitness.score)]
        setorder(unique.results, -fitness.score)
        final.result <- unique.results[seq_len(n.top.chroms), ]

        # print(paste('Algorithm terminated after', last.generation, 'generations.'))
        res.list <- list(top.chromosome.results = final.result, n.generations = last.generation)
        return(res.list)

    }
}

