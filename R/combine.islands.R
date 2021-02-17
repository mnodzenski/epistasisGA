#' A function to combine GADGETS results for individual islands into a single dataset.
#'
#' This function combines GADGETS results for individual islands into a single dataset.
#'
#' @param results.dir The directory in which individual island results from \code{run.gadgets} are saved.
#' @param annotation.data A data frame containing columns 'RSID', 'REF' and 'ALT'. Column 'RSID' gives the
#' RSIDs for the input SNPs, with the rows ordered such that the first RSID entry corresponds to the first SNP
#' column in the data passed to function \code{preprocess.genetic.data}, the second RSID corresponds to the second SNP column, etc.
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
#' library(Matrix)
#' block.ld.mat <- as.matrix(bdiag(list(matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25))))
#'
#' pp.list <- preprocess.genetic.data(case[, 1:10], father.genetic.data = dad[ , 1:10],
#'                                mother.genetic.data = mom[ , 1:10],
#'                                block.ld.mat = block.ld.mat[ , 1:10])
#'
#' run.gadgets(pp.list, n.chromosomes = 4, chromosome.size = 3, results.dir = 'tmp',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'
#' combined.res <- combine.islands('tmp', snp.annotations[ 1:10, ], pp.list)
#'
#' unlink('tmp', recursive = TRUE)
#' unlink('tmp_reg', recursive = TRUE)
#'
#' @importFrom data.table rbindlist setkey setorder `:=` setDT
#' @export

combine.islands <- function(results.dir, annotation.data, preprocessed.list) {

    # list all islands in the results data
    island.names <- list.files(results.dir, pattern = "cluster", full.names = TRUE)

    # note if we've already run this function
    out.file.name <- "combined.island.unique.chromosome.results.rds"
    out.file <- file.path(dirname(island.names[[1]]), out.file.name)
    if (file.exists(out.file)){

        message("combine.islands has already been run for this directory")

    }

    # grab the vector indicating which SNPs had coding flipped
    flipped.vec <- !preprocessed.list$minor.allele.vec

    # stop if the annotation data is not formatted correctly
    if (any(! c("RSID", "REF", "ALT") %in% colnames(annotation.data))){

        stop("annotation.data must contain columns RSID, REF, and ALT.")
    }

    if (nrow(annotation.data) != length(flipped.vec)){

        stop("annotation.data does not contain the same number of SNPs as the input data")
    }

    #flip the ref and alt alleles where we flipped the coding
    ref.current <- annotation.data$REF
    alt.current <- annotation.data$ALT
    annotation.data$REF[flipped.vec] <- alt.current[flipped.vec]
    annotation.data$ALT[flipped.vec] <- ref.current[flipped.vec]

    # then combine into a single data frame
    island.list <- lapply(island.names, function(island.file) {

        island <- gsub(".rds", "", basename(island.file))
        island.data <- readRDS(island.file)
        n.generations <- island.data$n.generations
        chrom.results <- island.data$top.chromosome.results
        chromosome.size <- sum(grepl("snp", colnames(chrom.results)))/3
        chrom.results[, `:=`(island, rep(island, nrow(chrom.results)))]
        chrom.results[, `:=`(n.generations, rep(n.generations, nrow(chrom.results)))]
        return(chrom.results)

    })

    # all results
    combined.result <- rbindlist(island.list)
    setorder(combined.result, -fitness.score)
    chromosome.size <- sum(grepl("snp", colnames(combined.result)))/3

    # subset to unique results
    combined.result[, `:=`(chromosome, paste(.SD, collapse = ".")), by = seq_len(nrow(combined.result)),
                    .SDcols = seq_len(chromosome.size)]
    unique.result <- combined.result[!duplicated(combined.result$chromosome), ]
    n.islands.found <- combined.result[, list(n.islands.found = length(fitness.score)), by = chromosome]
    setkey(unique.result, chromosome)
    setkey(n.islands.found, chromosome)
    unique.result <- unique.result[n.islands.found]
    unique.result[, `:=`(c("island", "n.generations"), NULL)]
    setorder(unique.result, -fitness.score)

    ## add in annotations for the SNPs and risk alleles
    risk.sign.cols <- seq_len(chromosome.size) + chromosome.size
    allele.copy.cols <- seq_len(chromosome.size) + 2*chromosome.size

    # starting with the rsids
    choose.these <- seq_len(chromosome.size)
    snp.cols <- unique.result[ , ..choose.these]
    snp.numbers <- unlist(snp.cols)
    rsids <- annotation.data$RSID
    rsid.dt <- data.table(matrix(rsids[snp.numbers], ncol = chromosome.size,
                                 byrow = FALSE))
    colnames(rsid.dt) <- paste(colnames(snp.cols), "rsid", sep = ".")

    #now the risk allele
    diff.cols <- unique.result[ , ..risk.sign.cols]
    diff.vecs <- unlist(diff.cols)
    risk.alleles <- rep(NA, length(diff.vecs))
    alt.alleles <- annotation.data$ALT
    ref.alleles <- annotation.data$REF
    risk.alleles[diff.vecs >= 0 ] <- alt.alleles[snp.numbers[diff.vecs >= 0]]
    risk.alleles[diff.vecs < 0 ] <- ref.alleles[snp.numbers[diff.vecs < 0]]
    risk.allele.dt <- data.table(matrix(risk.alleles, ncol = chromosome.size,
                                 byrow = FALSE))
    colnames(risk.allele.dt) <- gsub("diff.vec", "risk.allele", colnames(diff.cols))


    # put the full result together
    final.result <- cbind(snp.cols, rsid.dt, risk.allele.dt, unique.result[ , -(1:chromosome.size)])

    #write to file
    saveRDS(unique.result, file = out.file)

    # only unique chromosomes
    unique.file.name <- "combined.island.unique.chromosome.results.rds"
    unique.file <- file.path(dirname(island.names[[1]]), unique.file.name)
    saveRDS(final.result, file = unique.file)

    return(final.result)

}
