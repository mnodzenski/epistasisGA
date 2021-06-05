#' \code{epistasisGA} package
#'
#' A package implementing the GADGETS method to detect multi-SNP effects in case-parent triad or affected/unaffected sibling studies.
#'
#' @docType package
#' @name epistasisGAGE
NULL

## quiets concerns of R CMD check re visible bindings of global vars created via data table
if (getRversion() >= "2.15.1") {

    utils::globalVariables(c("chromosome", ".SD", "fitness.score", "raw.fitness.score",
        "Var1", "Var2", "name", "cluster", "job.id", "V1", "V2", "edge.score",
        "..choose.these", "..these.cols", "..rsid.cols", "SNP1", "SNP2", "SNP1.rsid",
        "SNP2.rsid", "NA_INTEGER", "..allele.copy.cols", "..risk.sign.cols","graphical.score",
        "data.type", "raw.score", "SNP", "pair.score", "rsid", "snp.score"))

}
