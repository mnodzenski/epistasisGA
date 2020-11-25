#' A function to plot a network of SNPs with potential multi-SNP effects.
#'
#' This function plots a network of SNPs with potential multi-SNP effects.
#'
#' @param edge.dt The data.table returned by function \code{compute.edge.scores}.
#' @param node.shape The desired node shape. See \code{names(igraph:::.igraph.shapes)} for available shapes.
#' @param repulse.rad A scalar affecting the graph shape. Decrease to reduce overlapping nodes,
#'  increase to move nodes closer together.
#' @param node.size A scalar affecting the size of the graph nodes. Increase to increase size.
#' @param graph.area A scalar affecting the size of the graph area. Increase to increase graph area.
#' @param vertex.label.cex A scalar controlling the size of the vertex label. Increase to increase size.
#' @param edge.width.cex A scalar controlling the width of the graph edges. Increase to make edges wider.
#' @param plot A logical indicating whether the network should be plotted. If set to false, this function will return an igraph object to be used for manual plotting.
#' @param edge.color.ramp A character vector of colors. The coloring of the network edges will be shown on a gradient, with the lower scoring edge weights
#' closer to the first color specified in \code{edge.color.ramp}, and higher scoring weights closer to the last color specified. By default, the
#' low scoring edges are whiter, and high scoring edges are redder.
#' @param node.color.ramp A character vector of colors. The coloring of the network nodes will be shown on a gradient, with the lower scoring nodes
#' closer to the first color specified in \code{node.color.ramp}, and higher scoring nodes closer to the last color specified. By default, the low
#' scoring nodes are whiter, and high scoring edges are greener.
#' @param ... Additional arguments to be passed to \code{plot.igraph}.
#' @return An igraph object, if \code{plot} is set to FALSE.
#'@examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' data(snp.annotations)
#' library(Matrix)
#' set.seed(1400)
#' block.ld.mat <- as.matrix(bdiag(list(matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25))))
#'
#' #preprocess data
#' target.snps <- c(1:3, 30:32, 60:62, 85)
#' pp.list <- preprocess.genetic.data(case[, target.snps], father.genetic.data = dad[ , target.snps],
#'                                mother.genetic.data = mom[ , target.snps],
#'                                block.ld.mat = block.ld.mat[target.snps , target.snps])
#' ## run GA for observed data
#'
#' #observed data chromosome size 2
#' run.ga(pp.list, n.chromosomes = 5, chromosome.size = 2, results.dir = 'tmp_2',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  combined.res2 <- combine.islands('tmp_2', snp.annotations[ target.snps, ], pp.list)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#'  #observed data chromosome size 3
#'  run.ga(pp.list, n.chromosomes = 5, chromosome.size = 3, results.dir = 'tmp_3',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  combined.res3 <- combine.islands('tmp_3', snp.annotations[ target.snps, ], pp.list)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#'  ## create list of results
#'  final.results <- list(combined.res2$unique.results[1:3, ], combined.res3$unique.results[1:3, ])
#'
#'  ## compute edge scores
#'  set.seed(20)
#'  edge.dt <- compute.edge.scores(final.results, pp.list, 3)
#'
#' ## plot
#' set.seed(10)
#' network.plot(edge.dt)
#'
#'  lapply(c('tmp_2', 'tmp_3'), unlink, recursive = TRUE)
#'
#' @import igraph
#' @importFrom qgraph qgraph.layout.fruchtermanreingold
#' @importFrom grDevices adjustcolor colorRampPalette
#' @importFrom data.table melt
#' @importFrom fields image.plot
#' @export

network.plot <- function(edge.dt, node.shape = "circle", repulse.rad = 1000,
    node.size = 25, graph.area = 100, vertex.label.cex = 0.5, edge.width.cex = 1, plot = TRUE,
    edge.color.ramp = c("white", "grey", "red"), node.color.ramp = c("yellow", "orange", "red"), ...) {

    #subset to target cols
    edge.dt <- edge.dt[ , c(3, 4, 5)]

    #compute node scores
    edge.dt.long <- melt(edge.dt, 3, c(1, 2), value.name = 'name')
    if (score.type == "max"){

        node.dt <- edge.dt.long[ , list(size = max(edge.score)), by = 'name']

    } else if (score.type == "sum"){

        node.dt <- edge.dt.long[ , list(size = sum(edge.score)), by = 'name']

    } else if (score.type == "logsum"){

        edge.dt.long[ , edge.score := exp(edge.score) - 1]
        node.dt <- edge.dt.long[ , list(size = log(1 + sum(edge.score))), by = 'name']

    }

    # convert to data.frames and scale the edge and node scores
    edge.df <- as.data.frame(edge.dt)
    edge.widths <- edge.df$edge.score/max(edge.df$edge.score)
    node.df <- as.data.frame(node.dt)
    node.df$size <- node.size*(node.df$size/max(node.df$size))

    # prepare for plotting
    colnames(edge.df)[1:2] <- c("from", "to")
    network <- graph.data.frame(edge.df[, 1:2], directed = FALSE, vertices = node.df)
    E(network)$weight <- edge.df$edge.score
    E(network)$width <- edge.width.cex*edge.widths
    color_fun <- colorRampPalette(edge.color.ramp)
    required.colors <- as.integer(as.factor(E(network)$weight))
    colors <- color_fun(length(unique(required.colors)))
    edge.colors <- sapply(seq_len(length(edge.widths)), function(x) adjustcolor(colors[required.colors][x],
        alpha.f = edge.widths[x]))
    E(network)$color <- edge.colors

    color_fun <- colorRampPalette(node.color.ramp)
    required.colors <- as.integer(as.factor(V(network)$size))
    colors <- color_fun(length(unique(required.colors)))
    V(network)$color <- colors[required.colors]
    V(network)$shape <- node.shape
    V(network)$label.cex <- vertex.label.cex*node.df$size/node.size

    # if desired, plot
    if (plot) {

        net.edges <- get.edgelist(network, names = FALSE)
        coords <- qgraph.layout.fruchtermanreingold(net.edges, vcount = vcount(network), repulse.rad = repulse.rad *
            vcount(network), area = graph.area * (vcount(network)^2))
        plot(network, layout = coords, asp = 0, ...)
        if (length(unique(edge.colors)) > 1){

            image.plot(legend.only = TRUE, zlim = range(V(network)$size), col = color_fun(500),
                       legend.lab = "SNP Score", legend.cex = 1.5, legend.line = 2.5)

        }

    # otherwise, return igraph object
    } else {

        return(network)

    }

}
