#' A function to combine GADGETS results for individual islands into a single dataset.
#'
#' This function combines GADGETS results for individual islands into a single dataset.
#'
#' @param results.dir The directory in which individual island results from \code{run.gadgets} are saved.
#' @param annotation.data A data frame containing columns 'RSID', 'REF' and 'ALT'. Column 'RSID' gives the
#' RSIDs for the input SNPs, with the rows ordered such that the first RSID entry corresponds to the first SNP
#' column in the data passed to function \code{preprocess.genetic.data}, the second RSID corresponds to the second SNP column, etc.
#' @param preprocessed.list The initial list produced by function \code{preprocess.genetic.data}.
#' @param n.top.chroms.per.island The number of top chromosomes per island to save in the final combined list. Defaults to the
#' top 10.
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
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'
#' combined.res <- combine.islands('tmp', snp.annotations[ 1:10, ], pp.list, 1)
#'
#' unlink('tmp', recursive = TRUE)
#' unlink('tmp_reg', recursive = TRUE)
#'
#' @importFrom data.table rbindlist setkey setorder `:=` setDT
#' @export

combine.islands <- function(results.dir, annotation.data, preprocessed.list, n.top.chroms.per.island = 1) {

    # list all islands in the results data
    island.names <- list.files(results.dir, pattern = "cluster", full.names = TRUE)

    # note if we've already run this function
    out.file.name <- "combined.island.unique.chromosome.results.rds"
    out.file <- file.path(dirname(island.names[1]), out.file.name)
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
        if (nrow(island.data$top.chromosome.results) < n.top.chroms.per.island ){

            stop("n.top.chroms.per.island must be <= the total number of chromosomes")

        }
        chrom.results <- island.data$top.chromosome.results

        # subset to unique results
        GxE <- "exposure" %in% colnames(chrom.results)
        chromosome.size <- sum(grepl("snp[0-9]$", colnames(chrom.results)))
        chrom.results[, `:=`(chromosome, paste(.SD, collapse = ".")), by = seq_len(nrow(chrom.results)),
                        .SDcols = seq_len(chromosome.size)]
        chrom.results <- chrom.results[!duplicated(chrom.results), ]

        #take top scorers
        chrom.results <- chrom.results[seq_len(n.top.chroms.per.island), ]
        return(chrom.results)

    })
    combined.result <- rbindlist(island.list)
    setorder(combined.result, -fitness.score)
    GxE <- any(grepl("exposure", colnames(combined.result)))
    chromosome.size <- sum(grepl("snp[0-9]$", colnames(combined.result)))

    # subset to unique results
    unique.result <- combined.result[!duplicated(combined.result$chromosome), ]
    n.islands.found <- combined.result[, list(n.islands.found = length(fitness.score)), by = chromosome]
    setkey(unique.result, chromosome)
    setkey(n.islands.found, chromosome)
    unique.result <- unique.result[n.islands.found]
    setorder(unique.result, -fitness.score)

    ## add in annotations for the SNPs and risk alleles

    # starting with the rsids
    choose.these <- seq_len(chromosome.size)
    snp.cols <- unique.result[ , ..choose.these]
    snp.numbers <- unlist(snp.cols)
    rsids <- annotation.data$RSID
    rsid.dt <- data.table(matrix(rsids[snp.numbers], ncol = chromosome.size,
                                 byrow = FALSE))
    colnames(rsid.dt) <- paste(colnames(snp.cols), "rsid", sep = ".")

    #now the risk allele
    alt.alleles <- annotation.data$ALT
    ref.alleles <- annotation.data$REF

    if (GxE){

        dif.vec.cols <- grep("^exposure.*diff.vec$", colnames(unique.result))
        allele.copy.cols <- grep("^exposure.*allele.copies$", colnames(unique.result))

        # loop over exposures and put together exposure specific results
        exposure.end.cols <- seq(chromosome.size, length(dif.vec.cols), chromosome.size)
        risk.allele.dt.list <- lapply(exposure.end.cols, function(end.col.idx){

            start.col.idx <- end.col.idx - chromosome.size + 1
            target.cols <- seq(dif.vec.cols[start.col.idx], dif.vec.cols[end.col.idx], 1)
            diff.cols <- unique.result[ , ..target.cols]
            diff.colnames <- colnames(diff.cols)
            diff.vecs <- unlist(diff.cols)
            risk.alleles <- rep(NA, length(diff.vecs))
            dv.gte0 <- !is.na(diff.vecs) & diff.vecs >= 0
            risk.alleles[dv.gte0] <- alt.alleles[snp.numbers[dv.gte0]]
            dv.lt0 <- !is.na(diff.vecs) & diff.vecs < 0
            risk.alleles[dv.lt0] <- ref.alleles[snp.numbers[dv.lt0]]
            risk.allele.dt <- data.table(matrix(risk.alleles, ncol = chromosome.size,
                                                byrow = FALSE))
            exposure.risk.allele.cols <- gsub("diff.vec", "risk.allele", diff.colnames)
            colnames(risk.allele.dt) <- exposure.risk.allele.cols
            exposure <- gsub(".snp1.diff.vec", "", diff.colnames[1])
            exposure.rank.col <- paste0(exposure, ".rank")
            exposure.risk.n.risk.allele.cols <- gsub("diff.vec", "allele.copies", diff.colnames)
            exposure.n.cases.risk.geno.col <- paste0(exposure, ".n.cases.risk.geno")
            exposure.n.comps.risk.geno.col <- paste0(exposure, ".n.comps.risk.geno")
            orig.target.cols <- c(exposure.rank.col, diff.colnames, exposure.risk.n.risk.allele.cols,
                                  exposure.n.cases.risk.geno.col, exposure.n.comps.risk.geno.col)
            exposure.original.dt <- unique.result[ , ..orig.target.cols]
            combined.dt <- cbind(exposure.original.dt, risk.allele.dt)
            new.target.cols <- c(exposure.rank.col, diff.colnames, exposure.risk.allele.cols,
                                 exposure.risk.n.risk.allele.cols, exposure.n.cases.risk.geno.col,
                                 exposure.n.comps.risk.geno.col)
            combined.dt <- combined.dt[ , ..new.target.cols]
            return(combined.dt)

        })

        # final results table
        overall.dif.vecs.and.fitness.cols <- c(paste0("snp", seq_len(chromosome.size), ".diff.vec"),
                                               "fitness.score")
        overall.dif.vecs.and.fitness <- unique.result[ , ..overall.dif.vecs.and.fitness.cols]
        exposure.specific.dt <- do.call(cbind, risk.allele.dt.list)
        n.islands <- unique.result[ , "n.islands.found"]
        final.result <- cbind(snp.cols, rsid.dt, overall.dif.vecs.and.fitness,
                              exposure.specific.dt, n.islands)

    } else {

        risk.sign.cols <- seq_len(chromosome.size) + chromosome.size
        diff.cols <- unique.result[ , ..risk.sign.cols]
        diff.vecs <- unlist(diff.cols)
        risk.alleles <- rep(NA, length(diff.vecs))
        risk.alleles[diff.vecs >= 0 ] <- alt.alleles[snp.numbers[diff.vecs >= 0]]
        risk.alleles[diff.vecs < 0 ] <- ref.alleles[snp.numbers[diff.vecs < 0]]
        risk.allele.dt <- data.table(matrix(risk.alleles, ncol = chromosome.size,
                                            byrow = FALSE))
        colnames(risk.allele.dt) <- gsub("diff.vec", "risk.allele", colnames(diff.cols))
        final.result <- cbind(snp.cols, rsid.dt, risk.allele.dt, unique.result[ , -(1:chromosome.size)])

    }

    # save
    saveRDS(final.result, file = out.file)

    return(final.result)

}

