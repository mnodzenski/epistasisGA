#' A function to combine results for individual islands into a single dataset.
#'
#' This function combines results for individual islands into a single dataset.
#'
#' @param results.dir The directory in which individual island results from \code{run.ga} are saved.
#' @return A data.table combing results from all islands stored in \code{results.dir}. This will also be written to \code{results.dir}
#' @importFrom data.table rbindlist
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

    island.data <- readRDS(island.file)
    top.chrom.results <- island.data$top.chromosome.results
    return(top.chrom.results)

  })
  combined.result <- rbindlist(island.list)
  setorder(combined.result, -fitness.score)
  out.file.name <- "combined.island.results.rds"
  out.file <- file.path(dirname(island.names[[1]]), out.file.name)
  saveRDS(combined.result, file = out.file)
  return(combined.result)

}
