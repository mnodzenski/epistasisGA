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
        chromosome.size <- sum(grepl("snp", colnames(chrom.results)))/3
        chrom.results[, `:=`(chromosome, paste(.SD, collapse = ".")), by = seq_len(nrow(chrom.results)),
                        .SDcols = seq_len(chromosome.size)]
        chrom.results <- chrom.results[!duplicated(chrom.results), ]

        #take top scorers
        chrom.results <- chrom.results[seq_len(n.top.chroms.per.island), ]
        return(chrom.results)

    })
    combined.result <- rbindlist(island.list)
    setorder(combined.result, -fitness.score)
    chromosome.size <- sum(grepl("snp", colnames(combined.result)))/3

    # subset to unique results
    unique.result <- combined.result[!duplicated(combined.result$chromosome), ]
    n.islands.found <- combined.result[, list(n.islands.found = length(fitness.score)), by = chromosome]
    setkey(unique.result, chromosome)
    setkey(n.islands.found, chromosome)
    unique.result <- unique.result[n.islands.found]
    setorder(unique.result, -fitness.score)

    ## add in annotations for the SNPs and risk alleles
<<<<<<< HEAD
    GxE <- "high.risk.exposure" %in% colnames(combined.result)
    if (GxE){

        chromosome.size <- sum(grepl("snp", colnames(combined.result)))/4
        risk.sign.cols <- seq_len(chromosome.size) + 2*chromosome.size
        allele.copy.cols <- seq_len(chromosome.size) + 3*chromosome.size

    } else {

        chromosome.size <- sum(grepl("snp", colnames(combined.result)))/3
        risk.sign.cols <- seq_len(chromosome.size) + chromosome.size
        allele.copy.cols <- seq_len(chromosome.size) + 2*chromosome.size

    }

    # starting with the rsids
    snp.col.positions <- seq_len(chromosome.size)
    snp.cols <- combined.result[ , ..snp.col.positions]
=======
    risk.sign.cols <- seq_len(chromosome.size) + chromosome.size
    allele.copy.cols <- seq_len(chromosome.size) + 2*chromosome.size

    # starting with the rsids
    choose.these <- seq_len(chromosome.size)
    snp.cols <- unique.result[ , ..choose.these]
>>>>>>> 8f61fd771ef587683168363a3c2e24e1dd581de1
    snp.numbers <- unlist(snp.cols)
    rsids <- annotation.data$RSID
    rsid.dt <- data.table(matrix(rsids[snp.numbers], ncol = chromosome.size,
                                 byrow = FALSE))
    colnames(rsid.dt) <- paste(colnames(snp.cols), "rsid", sep = ".")

    #now the risk allele
<<<<<<< HEAD
    diff.cols <- combined.result[ , ..risk.sign.cols]
=======
    diff.cols <- unique.result[ , ..risk.sign.cols]
>>>>>>> 8f61fd771ef587683168363a3c2e24e1dd581de1
    diff.vecs <- unlist(diff.cols)
    risk.alleles <- rep(NA, length(diff.vecs))
    alt.alleles <- annotation.data$ALT
    ref.alleles <- annotation.data$REF
    risk.alleles[diff.vecs >= 0 ] <- alt.alleles[snp.numbers[diff.vecs >= 0]]
    risk.alleles[diff.vecs < 0 ] <- ref.alleles[snp.numbers[diff.vecs < 0]]
    risk.allele.dt <- data.table(matrix(risk.alleles, ncol = chromosome.size,
                                 byrow = FALSE))
    if (GxE){

        colnames(risk.allele.dt) <- gsub("risk.sign", "risk.allele", colnames(diff.cols))

    } else {

        colnames(risk.allele.dt) <- gsub("diff.vec", "risk.allele", colnames(diff.cols))

    }

    ## count the number of cases and complements with the risk genotype
    original.col.numbers <- preprocessed.list$original.col.numbers
    case <- preprocessed.list$case.genetic.data
    comp <- preprocessed.list$complement.genetic.data

    n.case.comp.risk.geno.list <- lapply(seq_len(nrow(combined.result)), function(x){

        if (GxE){

            high.risk.exposure <- combined.result$high.risk.exposure[x]
            low.risk.exposure <- combined.result$low.risk.exposure[x]
            case.list <- lapply(c(high.risk.exposure, low.risk.exposure), function(exp.level){

                case[exposure == exp.level, ]

            })
            comp.list <- lapply(c(high.risk.exposure, low.risk.exposure), function(exp.level){

                comp[exposure == exp.level, ]

            })

        } else {

            case.list <- list(case)
            comp.list <- list(comp)

        }

        orig.chrom <- as.vector(t(snp.cols[x, ]))
        chrom <- which(original.col.numbers %in% orig.chrom)
        n.risk.alleles <- as.vector(t(combined.result[x, ..allele.copy.cols]))
        risk.signs <- sign(as.vector(t(diff.cols[x, ])))

        # determine the risk genotypes
        risk.geno <- ifelse(risk.signs >= 0 & n.risk.alleles == "2", 2,
                            ifelse(risk.signs < 0 & n.risk.alleles == "2", 0, 1))
        pos.cols <- risk.signs >= 0
        neg.cols <- risk.signs < 0

        # pick out the chromosome in the preprocessed list and the risk alleles
        unlist(lapply(seq_along(case.list), function(y){

            case <- case.list[[y]]
            comp <- comp.list[[y]]
            n <- nrow(case)

            # determine the number of cases and complements with the risk genotype
            risk.geno.mat <- matrix(rep(risk.geno, nrow(case)), nrow = nrow(case), byrow = TRUE)
            case.risk.geno <- matrix(NA, nrow = nrow(case), ncol = length(chrom))
            comp.risk.geno <- matrix(NA, nrow = nrow(case), ncol = length(chrom))

            if (any(pos.cols)){

                case.risk.geno[ , pos.cols] <- case[ , chrom[pos.cols]] >= risk.geno.mat[ , pos.cols]
                comp.risk.geno[ , pos.cols] <- comp[ , chrom[pos.cols]] >= risk.geno.mat[ , pos.cols]

            }
            if (any(neg.cols)){

                case.risk.geno[ , neg.cols] <- case[ , chrom[neg.cols]] >= risk.geno.mat[ , neg.cols]
                comp.risk.geno[ , neg.cols] <- comp[ , chrom[neg.cols]] >= risk.geno.mat[ , neg.cols]

            }
            n.case.full.risk.path <- sum(rowSums(case.risk.geno) == length(chrom))
            n.comp.full.risk.path <- sum(rowSums(comp.risk.geno) == length(chrom))
            if (GxE){

                return(c(n, n.case.full.risk.path, n.comp.full.risk.path))

            } else {

                return(c(n.case.full.risk.path, n.comp.full.risk.path))

            }

        }))

    })

    n.case.comp.risk.geno.dt <- t(setDT(n.case.comp.risk.geno.list))
    if (GxE){

        colnames(n.case.comp.risk.geno.dt) <- c("n.high.risk.exp", "n.cases.risk.geno.high.risk.exp", "n.comps.risk.geno.high.risk.exp",
                                                "n.low.risk.exp", "n.cases.risk.geno.low.risk.exp", "n.comps.risk.geno.low.risk.exp")

    } else {

        colnames(n.case.comp.risk.geno.dt) <- c("n.cases.risk.geno", "n.comps.risk.geno")
    }



    # put the full result together
<<<<<<< HEAD
    combined.result <- cbind(cbind(cbind(snp.cols, rsid.dt, risk.allele.dt), combined.result[ , -(1:chromosome.size)]),
                             n.case.comp.risk.geno.dt)
    combined.result[, `:=`(chromosome, paste(.SD, collapse = ".")), by = seq_len(nrow(combined.result)),
                  .SDcols = seq_len(chromosome.size)]

    #write to file
    saveRDS(combined.result, file = out.file)
=======
    final.result <- cbind(snp.cols, rsid.dt, risk.allele.dt, unique.result[ , -(1:chromosome.size)])
>>>>>>> 8f61fd771ef587683168363a3c2e24e1dd581de1

    # save
    saveRDS(final.result, file = out.file)

    return(final.result)

}
