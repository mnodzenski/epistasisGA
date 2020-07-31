#' A function to plot a network of SNPs with potential multi-SNP effects.
#'
#' This function plots a network of SNPs with potential multi-SNP effects.
#'
#' @param results.df A subset of the \code{unique.results} data frame of results from \code{combine.islands} after running \code{run.ga}.If not
#' specified, \code{edge.dt} must be specified.
#' @param edge.dt The data.table returned by function \code{network.threshold}. If not specified, \code{results.df} must instead be specified.
#' @param node.shape The desired node shape. See \code{names(igraph:::.igraph.shapes)} for available shapes.
#' @param score.type A character string specifying the method for weighting network edges and nodes, with options
#' 'max', 'sum', or 'logsum'. The default is 'max', but 'logsum' may also be particularly useful.
#'  Note that "logsum" is actually the log of one plus the sum of the fitness scores to avoid nodes or edges having negative
#'  weights.
#' @param repulse.rad A scalar affecting the graph shape. Decrease to reduce overlapping nodes, increase to move nodes closer together.
#' @param node.size A scalar affecting the size of the graph nodes. Increase to increase size.
#' @param graph.area A scalar affecting the size of the graph area. Increase to increase graph area.
#' @param vertex.label.cex A scalar controlling the size of the vertex label. Increase to increase size.
#' @param plot A logical indicating whether the network should be plotted. If set to false, this function will return an igraph object to be used for manual plotting.
#' @return An igraph object, if \code{plot} is set to FALSE.
#'@examples
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
#' run.ga(pp.list, n.chromosomes = 4, chromosome.size = 3, results.dir = 'tmp',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'
#' combined.res <- combine.islands('tmp', snp.annotations[ 1:10, ], pp.list)
#'
#' set.seed(10)
#' network.plot(combined.res$unique.results)
#'
#' unlink('tmp', recursive = TRUE)
#' unlink('tmp_reg', recursive = TRUE)
#'
#' @import igraph
#' @importFrom qgraph qgraph.layout.fruchtermanreingold
#' @importFrom grDevices adjustcolor colorRampPalette
#' @importFrom data.table melt
#' @export

network.plot <- function(results.df = NULL, edge.dt = NULL, node.shape = "crectangle", score.type = "max", repulse.rad = 1000,
    node.size = 25, graph.area = 100, vertex.label.cex = 0.5, plot = TRUE) {

    # if not inputting an edge.df, compute it
    if (is.null(edge.dt)){

        edge.dt <- compute.edge.scores(results.df, score.type = score.type)

    }

    #compute node scores
    edge.dt.long <- melt(edge.dt, 3, 1:2, value.name = 'name')
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
    E(network)$width <- edge.widths
    color_fun <- colorRampPalette(c("white", "grey", "red"))
    required.colors <- as.integer(as.factor(E(network)$weight))
    colors <- color_fun(length(unique(required.colors)))
    edge.colors <- sapply(seq_len(length(edge.widths)), function(x) adjustcolor(colors[required.colors][x],
        alpha.f = edge.widths[x]))
    E(network)$color <- edge.colors

    color_fun <- colorRampPalette(c("white", "grey", "green"))
    required.colors <- as.integer(as.factor(V(network)$size))
    colors <- color_fun(length(unique(required.colors)))
    V(network)$color <- colors[required.colors]
    V(network)$shape <- node.shape

    # if desired, plot
    if (plot) {

        net.edges <- get.edgelist(network, names = FALSE)
        coords <- qgraph.layout.fruchtermanreingold(net.edges, vcount = vcount(network), repulse.rad = repulse.rad *
            vcount(network), area = graph.area * (vcount(network)^2))
        plot(network, layout = coords, vertex.label.cex = vertex.label.cex, asp = 0)

    # otherwise, return igraph object
    } else {

        return(network)

    }

}
