#' A function to combine GADGETS results for individual islands into a single
#' dataset.
#'
#' This function combines GADGETS results for individual islands into a single
#' dataset.
#'
#' @param results.dir The directory in which individual island results from
#' \code{run.gadgets} are saved.
#' @param annotation.data A data frame containing columns 'RSID', 'REF' and
#' 'ALT'. Column 'RSID' gives the RSIDs for the input SNPs, with the rows
#' ordered such that the first RSID entry corresponds to the first SNP
#' column in the data passed to function \code{preprocess.genetic.data}, the
#' second RSID corresponds to the second SNP column, etc.
#' @param preprocessed.list The initial list produced by function
#' \code{preprocess.genetic.data}.
#' @param n.top.chroms.per.island The number of top chromosomes per island to
#' save in the final combined list. Defaults to the single top chromosome.
#' @return A data.table containing the results aggregated across islands. Note
#' these results be written to \code{results.dir} as
#' combined.island.unique.chromosome.results.rds'. See the package vignette for
#' more detailed descriptions of the content of each output column. Secondarily,
#' this will concatenate all individual island results files and store them
#' in a single file, called "all.island.results.concatenated.rds".
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' data(snp.annotations)
#'
#' pp.list <- preprocess.genetic.data(as.matrix(case[, 1:10]),
#'                                father.genetic.data = as.matrix(dad[ , 1:10]),
#'                                mother.genetic.data = as.matrix(mom[ , 1:10]),
#'                                ld.block.vec = c(10))
#'
#' run.gadgets(pp.list, n.chromosomes = 4, chromosome.size = 3,
#'        results.dir = 'tmp',
#'        cluster.type = 'interactive',
#'        registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1,
#'        n.migrations = 0)
#'
#' combined.res <- combine.islands('tmp', snp.annotations[ 1:10, ], pp.list)
#'
#' unlink("tmp", recursive = TRUE)
#' unlink("tmp_reg", recursive = TRUE)
#'
#' @importFrom data.table rbindlist setkey setorder `:=` setDT
#' @export

combine.islands <- function(results.dir, annotation.data, preprocessed.list,
                            n.top.chroms.per.island = 1) {

    # list all islands in the results data
    island.names <- list.files(results.dir, pattern = "cluster",
                               full.names = TRUE)

    # note if we've already run this function
    out.file.name <- "combined.island.unique.chromosome.results.rds"
    out.file <- file.path(dirname(island.names[1]), out.file.name)
    if (file.exists(out.file)) {
        message("combine.islands has already been run for this directory")
    }

    # stop if the annotation data is not formatted correctly
    if (any(!c("RSID", "REF", "ALT") %in% colnames(annotation.data))) {
        stop("annotation.data must contain columns RSID, REF, and ALT.")
    }

    n.candidate.snps <- ncol(preprocessed.list$case.genetic.data)
    if (nrow(annotation.data) != n.candidate.snps){

        stop("annotation.data does not contain the same number of SNPs as the input data")
    }

    # then combine into a single data frame
    island.list <- lapply(island.names, function(island.file) {
        island <- gsub(".rds", "", basename(island.file))
        island.data <- readRDS(island.file)
        n.generations <- island.data$n.generations
        if (nrow(island.data$top.chromosome.results) < n.top.chroms.per.island)
            {
            stop("n.top.chroms.per.island must be <= the total chromosomes")
        }
        chrom.results <- island.data$top.chromosome.results

        # subset to unique results
        chromosome.size <- sum(grepl("snp[0-9]$", colnames(chrom.results)))
        chrom.results[, `:=`(chromosome, paste(.SD, collapse = ".")),
                      by = seq_len(nrow(chrom.results)),
                        .SDcols = seq_len(chromosome.size)]

        # also saving the full results for combined file
        all.results <- chrom.results
        all.results$island <- island
        all.results$n.generations <- n.generations
        chrom.results <- chrom.results[!duplicated(chrom.results), ]

        # take top scorers
        chrom.results <- chrom.results[seq_len(n.top.chroms.per.island), ]
        return(list(chrom.results, all.results))

    })
    combined.result <- rbindlist(lapply(island.list, function(x) x[[1]]))
    setorder(combined.result, -fitness.score)

    all.island.res <- rbindlist(lapply(island.list, function(x) x[[2]]))
    
    # rename the exposure level betas based on the input data 
    if (preprocessed.list$E_GADGETS){
        
        in.exp.cols <- paste0(colnames(preprocessed.list$exposure.mat), 
                              "_p_disease_coef")
        these.cols <- grepl("risk.exp.beta", colnames(combined.result), 
                            fixed = TRUE)
        new.names <- c("intercept_p_disease_coef", in.exp.cols)
        colnames(combined.result)[these.cols] <- new.names        
        colnames(all.island.res)[these.cols] <- new.names
        
    }

    #remove all the individual island files
    lapply(island.names, unlink)

    # save concatenated file instead
    all.islands.out.file <- file.path(dirname(island.names[1]),
                                      "all.island.results.concatenated.rds")
    saveRDS(all.island.res, all.islands.out.file)
    chromosome.size <- sum(grepl("snp[0-9]$", colnames(combined.result)))

    # subset to unique results
    unique.result <- combined.result[!duplicated(combined.result$chromosome), ]
    n.islands.found <- combined.result[, list(n.islands.found =
                                                  length(fitness.score)),
                                       by = chromosome]
    setkey(unique.result, chromosome)
    setkey(n.islands.found, chromosome)
    unique.result <- unique.result[n.islands.found]
    setorder(unique.result, -fitness.score)

    ## add in annotations for the SNPs and risk alleles

    # starting with the rsids
    choose.these <- seq_len(chromosome.size)
    snp.cols <- unique.result[, ..choose.these]
    snp.numbers <- unlist(snp.cols)
    rsids <- annotation.data$RSID
    rsid.dt <- data.table(matrix(rsids[snp.numbers],
        ncol = chromosome.size,
        byrow = FALSE
    ))
    colnames(rsid.dt) <- paste(colnames(snp.cols), "rsid", sep = ".")

    #now the risk allele
    alt.alleles <- annotation.data$ALT
    ref.alleles <- annotation.data$REF
    risk.sign.cols <- seq_len(chromosome.size) + chromosome.size
    diff.cols <- unique.result[ , ..risk.sign.cols]
    diff.vecs <- unlist(diff.cols)
    risk.alleles <- rep(NA, length(diff.vecs))
    risk.alleles[diff.vecs >= 0 ] <- alt.alleles[
        snp.numbers[diff.vecs >= 0]]
    risk.alleles[diff.vecs < 0 ] <- ref.alleles[snp.numbers[diff.vecs < 0]]
    risk.allele.dt <- data.table(matrix(risk.alleles,
                                        ncol = chromosome.size,
                                        byrow = FALSE))
    colnames(risk.allele.dt) <- gsub("diff.vec", "risk.allele",
                                     colnames(diff.cols))
    not.these <- -seq_len(chromosome.size)
    final.result <- cbind(snp.cols, rsid.dt, risk.allele.dt,
                          unique.result[ , ..not.these])
    
    # save
    saveRDS(final.result, file = out.file)
    return(final.result)
}

