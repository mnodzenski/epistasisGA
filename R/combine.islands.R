#' A function to combine results for individual islands into a single dataset.
#'
#' This function combines results for individual islands into a single dataset.
#'
#' @param results.dir The directory in which individual island results from \code{run.ga} are saved.
#' @return A list of two elements, \code{all.results} and \code{unique.results}. \code{all.results} contains chromosome results across all islands, including chromsomes identified on more than one island. \code{unique.results} removes duplicate chromosomes found on more than one island. These two objects will also be written to \code{results.dir} as combined.island.results.rds and combined.island.unique.chromosome.results.rds.
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' library(Matrix)
#' chrom.mat <- as.matrix(bdiag(list(matrix(rep(TRUE, 2500^2), nrow = 2500),
#'                               matrix(rep(TRUE, 2500^2), nrow = 2500),
#'                               matrix(rep(TRUE, 2500^2), nrow = 2500),
#'                               matrix(rep(TRUE, 2500^2), nrow = 2500))))
#'
#' pp.list <- preprocess.genetic.data(case[, 1:10], father.genetic.data = dad[ , 1:10],
#'                                mother.genetic.data = mom[ , 1:10],
#'                                chrom.mat = chrom.mat[ , 1:10])
#'
#' run.ga(pp.list, n.chromosomes = 4, chromosome.size = 3, results.dir = "tmp",
#'        cluster.type = "interactive", registryargs = list(file.dir = "tmp_reg", seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'
#' combined.res <- combine.islands("tmp")
#'
#' unlink("tmp", recursive = TRUE)
#' unlink("tmp_reg", recursive = TRUE)
#'
#' @importFrom data.table rbindlist setkey setorder
#' @export

combine.islands <- function(results.dir){

  #list all islands in the results data
  island.names <- list.files(results.dir, full.names = TRUE)

  #stop if the islands have already been combined
  if (any(grepl("combined.island", island.names))){

    stop("Islands have already been combined")

  }

  #otherwise combine into a single data frame
  island.list <- lapply(island.names, function(island.file){

    island <- gsub(".rds", "", basename(island.file))
    island.data <- readRDS(island.file)
    n.generations <- island.data$n.generations
    chrom.results <- island.data$top.chromosome.results
    chromosome.size <- sum(grepl("snp", colnames(chrom.results)))/2
    chrom.results[ , island := rep(island, nrow(chrom.results))]
    chrom.results[ , n.generations := rep(n.generations, nrow(chrom.results))]
    chrom.results[ , chromosome := paste(.SD, collapse = "."), by = seq_len(nrow(chrom.results)), .SDcols = 1:chromosome.size]
    return(chrom.results)

  })

  #all results
  combined.result <- rbindlist(island.list)
  setorder(combined.result, -fitness.score)
  out.file.name <- "combined.island.results.rds"
  out.file <- file.path(dirname(island.names[[1]]), out.file.name)
  saveRDS(combined.result, file = out.file)

  #only unique chromosomes
  unique.result <- combined.result[!duplicated(combined.result$chromosome), ]
  n.islands.found <- combined.result[ , list(n.islands.found = length(fitness.score)), by = chromosome]
  setkey(unique.result, chromosome)
  setkey(n.islands.found, chromosome)
  unique.result <- unique.result[n.islands.found]
  unique.result[ , c("island", "n.generations") := NULL]
  setorder(unique.result, -fitness.score)
  unique.file.name <- "combined.island.unique.chromosome.results.rds"
  unique.file <- file.path(dirname(island.names[[1]]), unique.file.name)
  saveRDS(unique.result, file = unique.file)

  return(list(all.results = combined.result, unique.results = unique.result))

}
