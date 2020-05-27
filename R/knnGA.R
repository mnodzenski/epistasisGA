#' \code{knnGA} package
#'
#' A package implementing a genetic algorithm to detect multi-SNP effects in case-parent triad studies.
#'
#' @docType package
#' @name knnGA
NULL

## quiets concerns of R CMD check re visible bindings of global vars created via data table
if (getRversion() >= "2.15.1") {
    
    utils::globalVariables(c("chromosome", ".SD", "fitness.score", "raw.fitness.score", "min.elem", 
        "Var1", "Var2", "name", "cluster", "job.id"))
    
}
