#' A function to combine results for individual islands into a single dataset.
#'
#' This function combines results for individual islands into a single dataset. It can only be run after calling function \code{combine.islands}.
#'
#' @param results.dir The directory in which individual island results from \code{run.ga} are saved.
#' @param preprocessed.list The initial list produced by function \code{preprocess.genetic.data}.
#' @param n.different.snps.weight The number by which the number of different SNPs between a case and complement is multiplied in computing the family weights. Defaults to 2.
#' @param n.both.risk.geno.weight The number by which the number of SNPs with the risk genotype in both the case and complement is multiplied in computing the family weights. Defaults to 1.
#' @param weight.function.int An integer used to assign family weights. Specifically, we use \code{weight.function.int} in a  function that takes the weighted sum
#' of the number of different SNPs and SNPs both equal to one as an argument, denoted as x, and returns a family weight equal to \code{weight.function.int}^x. Defaults to 2.
#' @return A list of two elements. Note these two objects will also be written to \code{results.dir}
#' as 'epi.combined.island.results.rds' and 'epi.combined.island.unique.chromosome.results.rds'. Furthermore,
#' this function should not be called twice on the same directory.
#' \describe{
#'  \item{epi.all.results}{A dataset containing chromosome results across all islands,
#'  where top chromosomes that evolved on multiple distinct islands appear in multiple rows.}
#'  \item{epi.unique.results}{A condensed version of \code{all.results} with one row per distinct chromosome
#'  and an additional variable indicating the number of islands on which that chromosome evolved.}
#' }
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' data(snp.annotations)
#' library(Matrix)
#' chrom.mat <- as.matrix(bdiag(list(matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25))))
#'
#' pp.list <- preprocess.genetic.data(case[, 1:10], father.genetic.data = dad[ , 1:10],
#'                                mother.genetic.data = mom[ , 1:10],
#'                                chrom.mat = chrom.mat[ , 1:10])
#'
#' run.ga(pp.list, n.chromosomes = 4, chromosome.size = 3, results.dir = 'tmp',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'
#' combined.res <- combine.islands('tmp', snp.annotations[ 1:10, ], pp.list)
#' epi.score.res <- compute.epi.fitness.scores('tmp', pp.list)
#'
#' unlink('tmp', recursive = TRUE)
#' unlink('tmp_reg', recursive = TRUE)
#'
#' @importFrom data.table data.table rbindlist setkey setorder `:=`
#' @export

compute.epi.fitness.scores <- function(results.dir, preprocessed.list, n.different.snps.weight = 2,
                                       n.both.risk.geno.weight = 1,  weight.function.int = 2) {

    # list all islands in the results data
    island.names <- list.files(results.dir, full.names = TRUE)

    # make sure we have already combined.islands
    ci.file.name <- "combined.island.results.rds"
    ci.file <- file.path(dirname(island.names[[1]]), ci.file.name)
    if (!file.exists(ci.file)){

        stop("combine.islands has not yet been run for this directory")

    }

    # grab the previous unique results
    prev.res <- readRDS(file.path(dirname(island.names[[1]]), 'combined.island.unique.chromosome.results.rds'))

    # determine chromosome size
    chrom.size <- sum(grepl("snp", colnames(prev.res)))/5

    # set up weight lookup
    max.sum <- max(n.different.snps.weight, n.both.risk.geno.weight)*chrom.size
    weight.lookup <- vapply(seq_len(max.sum), function(x) weight.function.int^x, 1)

    # recompute fitness scores using epistasis weighting
    new.scores <- rbindlist(lapply(seq_len(nrow(prev.res)), function(res.row){

        # pick out the relevant information for each row of the results
        snps <- as.vector(t(prev.res[res.row, 1:chrom.size]))
        diff.vec.cols <- paste0("snp", 1:chrom.size, ".diff.vec")
        diff.vec <- as.vector(t(prev.res[res.row, ..diff.vec.cols]))
        allele.copy.cols <- paste0("snp", 1:chrom.size, ".allele.copies")
        allele.copies <- as.vector(t(prev.res[res.row, ..allele.copy.cols]))

        # map the original columns to the preprocessed list column
        target.cols <- match(snps, preprocessed.list$original.col.numbers)
        epi.res <- epi.chrom.fitness.score(preprocessed.list$case.genetic.data,
                                preprocessed.list$complement.genetic.data,
                                target.cols, diff.vec, allele.copies,
                                preprocessed.list$chrom.mat, weight.lookup)
        raw.epi.score <- epi.res$fitness.score
        epi.score <- min(abs(epi.res$sum.dif.vecs))*raw.epi.score
        warn.message <- epi.res$warn.message

        # return result
        res <- data.table(raw.epi.fitness.score = raw.epi.score, epi.fitness.score = epi.score,
                          warnings = warn.message)

    }))

    # append the new scores
    updated.res <- cbind(prev.res, new.scores)

    # sort by epi.fitness.score
    setorder(updated.res, -epi.fitness.score)

    # generate the non-unique version of the data
    full.updated.res <- updated.res[rep(seq_len(nrow(updated.res)), updated.res$n.islands.found)]

    #write to file
    saveRDS(updated.res,
            file = file.path(dirname(island.names[[1]]), 'epi.combined.island.unique.chromosome.results.rds'))
    saveRDS(full.updated.res,
            file = file.path(dirname(island.names[[1]]), 'combined.island.unique.chromosome.results.rds'))

    return(list(all.results = full.updated.res, unique.results = updated.res))

}
