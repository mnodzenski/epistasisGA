#' A function to compute edge scores for network plots of results.
#'
#' This function returns a data.table of edge weights for use in network plots of GA results.
#'
#' @param results.df A subset of the \code{unique.results} data frame of results from \code{combine.islands} after running \code{run.ga}.
#' @param score.type A character string specifying the method for scoring edges, with options
#' 'max', 'sum', or 'logsum'. The default is 'max', but 'logsum' may also be particularly useful.
#'  Note that "logsum" is actually the log of one plus the sum of the fitness scores to avoid nodes or edges having negative
#'  weights.
#' @return A data.table where the first two columns represent SNPs and the third column (edge.score)
#' is the edge score of a chromosome containing those SNPs.
#'
#'@examples
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
#'
#' set.seed(10)
#' edge.scores <- compute.edge.scores(combined.res$unique.results)
#'
#' unlink('tmp', recursive = TRUE)
#' unlink('tmp_reg', recursive = TRUE)
#'
#' @importFrom data.table data.table rbindlist setorder
#' @importFrom utils combn
#' @export

compute.edge.scores <- function(results.df, score.type = "max") {

    chrom.size <- sum(grepl("snp", colnames(results.df)))/5
    n.top.chroms <- nrow(results.df)
    all.edge.weights <- rbindlist(lapply(seq_len(n.top.chroms), function(res.row) {

        chrom.res <- results.df[res.row, ]
        fs <- chrom.res$fitness.score
        hs <- chrom.res$fitness.score*chrom.res$n.islands.found
        chrom <- as.vector(t(chrom.res[, 1:chrom.size]))
        chrom.pairs <- data.table(t(combn(chrom, 2)))
        chrom.pairs[ , `:=`(fitness.score = fs, h.score = hs)]
        return(chrom.pairs)

    }))

    # figure out the edge score based on score.type
    if (score.type == "max"){

        out.dt <- all.edge.weights[ , .(edge.score = max(fitness.score)), .(V1, V2)]
        setorder(out.dt, -edge.score)

    } else if (score.type == "sum"){

        out.dt <- all.edge.weights[ , .(edge.score = sum(h.score)), .(V1, V2)]
        setorder(out.dt, -edge.score)

    } else if (score.type == "logsum"){

        out.dt <- all.edge.weights[ , .(edge.score = log(1 + sum(h.score))), .(V1, V2)]
        setorder(out.dt, -edge.score)

    }

    return(out.dt)

}


