#' \code{epistasisGA} package
#'
#' A package implementing the GADGETS method to detect multi-SNP effects in case-parent triad studies.
#'
#' @docType package
#' @name epistasisGA
NULL

## quiets concerns of R CMD check re visible bindings of global vars created via data table
if (getRversion() >= "2.15.1") {

    utils::globalVariables(c("chromosome", ".SD", "fitness.score", "raw.fitness.score", "min.elem",
        "Var1", "Var2", "name", "cluster", "job.id", "V1", "V2", "edge.score", "h.score",
        "..choose.these"))

}
