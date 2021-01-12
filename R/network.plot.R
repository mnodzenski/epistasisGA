#' A function to plot a network of SNPs with potential multi-SNP effects.
#'
#' This function plots a network of SNPs with potential multi-SNP effects.
#'
#' @param edge.dt The data.table returned by function \code{compute.edge.scores}, or a subset of it. By default, the SNPs
#' will be labeled with their RSIDs, listed in columns 3 and 4. Users can create custom labels by changing the values in these
#' two columns.
#' @param preprocessed.list The initial list produced by function \code{preprocess.genetic.data}.
#' @param score.type A character string specifying the method for aggregating SNP-pair scores across chromosome sizes. Options are
#' 'max', 'sum', or 'logsum', defaulting to "logsum". For a given SNP-pair, it's graphical score will be the \code{score.type} of all
#' graphical scores of chromosomes containing that pair across chromosome sizes. Pair scores will be proportional to the sum of graphical scores
#' for either 'logsum' or 'sum', but 'logsum' may be useful in cases where there are multiple risk-sets, and one is found much more frequently.
#' Note that "logsum" is actually the log of one plus the sum of the SNP-pair scores to avoid nodes or edges having negative weights.
#' @param node.shape The desired node shape. See \code{names(igraph:::.igraph.shapes)} for available shapes. Defaults to circle.
#' @param repulse.rad A scalar affecting the graph shape. Decrease to reduce overlapping nodes,
#'  increase to move nodes closer together.
#' @param node.size A scalar affecting the size of the graph nodes. Increase to increase size.
#' @param graph.area A scalar affecting the size of the graph area. Increase to increase graph area.
#' @param vertex.label.cex A scalar controlling the size of the vertex label. Increase to increase size.
#' @param edge.width.cex A scalar controlling the width of the graph edges. Increase to make edges wider.
#' @param plot A logical indicating whether the network should be plotted. If set to false, this function will return an igraph object to be used for manual plotting.
#' @param edge.color.ramp A character vector of colors. The coloring of the network edges will be shown on a gradient, with the lower scoring edge weights
#' closer to the first color specified in \code{edge.color.ramp}, and higher scoring weights closer to the last color specified. By default, the
#' low scoring edges are light blue, and high scoring edges are dark blue.
#' @param node.color.ramp A character vector of colors. The coloring of the network nodes will be shown on a gradient, with the lower scoring nodes
#' closer to the first color specified in \code{node.color.ramp}, and higher scoring nodes closer to the last color specified. By default, the low
#' scoring nodes are whiter, and high scoring edges are redder.
#' @param plot.legend A boolean indicating whether a legend should be plotted. Defaults to TRUE.
#' @param high.ld.threshold A numeric value between 0 and 1, indicating the r^2 threshold in complements (or unaffected siblings)
#' above which a pair of SNPs in the same LD block (as specified in \code{preprocessed.list}) should be considered in high LD. Connections
#' between these high LD SNPs will be dashed instead of solid lines. Defaults to 0.25.
#' @param plot.margins A vector of length 4 passed to \code{par(mar = )}. Defaults to c(2, 1, 2, 1).
#' @param legend.title.cex A numeric value controlling the size of the legend titles. Defaults to 1.75. Increase
#' to increase font size, decrease to decrease font size.
#' @param legend.axis.cex A numeric value controlling the size of the legend axis labels. Defaults to 1.75. Increase
#' to increase font size, decrease to decrease font size.
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
#' run.gadgets(pp.list, n.chromosomes = 5, chromosome.size = 2, results.dir = 'tmp_2',
#'        cluster.type = 'interactive', registryargs = list(file.dir = 'tmp_reg', seed = 1500),
#'        generations = 2, n.islands = 2, island.cluster.size = 1, n.top.chroms = 3)
#'  combined.res2 <- combine.islands('tmp_2', snp.annotations[ target.snps, ], pp.list)
#'  unlink('tmp_reg', recursive = TRUE)
#'
#'  #observed data chromosome size 3
#'  run.gadgets(pp.list, n.chromosomes = 5, chromosome.size = 3, results.dir = 'tmp_3',
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
#'  edge.dt <- compute.pair.scores(final.results, pp.list, 3, pval.thresh = 1)
#'
#' ## plot
#' set.seed(10)
#' network.plot(edge.dt, pp.list)
#'
#'  lapply(c('tmp_2', 'tmp_3'), unlink, recursive = TRUE)
#'
#' @import igraph
#' @importFrom qgraph qgraph.layout.fruchtermanreingold
#' @importFrom grDevices adjustcolor colorRampPalette as.raster
#' @importFrom data.table melt
#' @importFrom stats cor
#' @importFrom graphics rasterImage axis layout par
#' @export

network.plot <- function(edge.dt, preprocessed.list, score.type = "logsum", node.shape = "circle",
                         repulse.rad = 1000, node.size = 25, graph.area = 100, vertex.label.cex = 0.5,
                         edge.width.cex = 1, plot = TRUE, edge.color.ramp = c("lightblue", "blue"),
                         node.color.ramp = c("white", "red"), plot.legend = TRUE,
                         high.ld.threshold = 0.1, plot.margins = c(2, 1, 2, 1), legend.title.cex = 1.75,
                         legend.axis.cex = 1.75, ...) {

    #compute r2 vals for snps in the same ld block, assign 0 otherwise
    original.col.numbers <- preprocessed.list$original.col.numbers
    r2.vals <- vapply(seq(1, nrow(edge.dt)), function(x){

        # pick out the snp pair in the preprocessed list
        snp.pair <- as.vector(t(edge.dt[x, c(1, 2)]))
        target.snps <- which(original.col.numbers %in% snp.pair)

        # check if snps are located in same ld block
        block.ld.mat <- preprocessed.list$block.ld.mat
        target.block.ld.mat <- block.ld.mat[target.snps, target.snps]
        same.ld.block <- target.block.ld.mat[2, 1]

        # if not on same ld block, compute r2
        if (!same.ld.block){

            return(0.0)

        } else {

            comp.genetic.data <- preprocessed.list$complement.genetic.data
            snp1 <- target.snps[1]
            snp2 <- target.snps[2]
            r2 <- cor(comp.genetic.data[ , snp1], comp.genetic.data[ , snp2])^2
            return(r2)

        }

    }, 1.0)


    #subset to target cols
    edge.label.dt <- edge.dt[ , c(3, 4, 5)]
    edge.dt <- edge.dt[ , c(1, 2, 5)]

    #compute node scores
    edge.dt.long <- melt(edge.dt, 3, c(1, 2), value.name = 'name')
    edge.label.dt.long <- melt(edge.label.dt, 3, c(1, 2), value.name = 'label')
    node.labels <- as.character(edge.label.dt.long$label)
    names(node.labels) <- as.character(edge.dt.long$name)

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
    node.size.raw <- node.df$size/max(node.df$size)
    node.df$size <- node.size*(node.df$size/max(node.df$size))

    # prepare for plotting
    colnames(edge.df)[1:2] <- c("from", "to")
    network <- graph.data.frame(edge.df[, 1:2], directed = FALSE, vertices = node.df)
    E(network)$weight <- edge.df$edge.score
    E(network)$width <- edge.width.cex*edge.widths
    color_fun_e <- colorRampPalette(edge.color.ramp)
    edge.required.colors <- as.integer(as.factor(E(network)$weight))
    raw.edge.colors <- color_fun_e(length(unique(edge.required.colors)))
    edge.colors <- sapply(seq_len(length(edge.widths)),
                          function(x) adjustcolor(raw.edge.colors[edge.required.colors][x],
                                                  alpha.f = edge.widths[x]))
    E(network)$color <- edge.colors

    color_fun_n <- colorRampPalette(node.color.ramp)
    node.required.colors <- as.integer(as.factor(V(network)$size))
    node.colors <- color_fun_n(length(unique(node.required.colors)))
    V(network)$color <- node.colors[node.required.colors]
    V(network)$shape <- node.shape
    V(network)$label.cex <- vertex.label.cex*node.df$size/node.size
    V(network)$label <- node.labels[V(network)$name]

    E(network)$lty <- ifelse(r2.vals >= high.ld.threshold, 3, 1)

    # if desired, plot
    if (plot) {

        net.edges <- get.edgelist(network, names = FALSE)
        coords <- qgraph.layout.fruchtermanreingold(net.edges, vcount = vcount(network), repulse.rad = repulse.rad *
            vcount(network), area = graph.area * (vcount(network)^2))

        if (length(unique(edge.colors)) > 1 & plot.legend){

            par(mar = plot.margins)
            layout(matrix(c(1, 1, 2, 3), ncol = 2, byrow = F), widths = c(3.5,0.5), heights = c(1,1))
            plot(network, layout = coords, asp = 0, ...)

            node_legend <- as.raster(matrix(rev(node.colors), ncol = 1))
            plot(c(0,2),c(0,1),type = 'n', axes = F, xlab = '', ylab = '', main = 'SNP-Score',
                 cex.main = legend.title.cex)
            rasterImage(node_legend, 0.75, 0, 1, 1)
            n.legend.labels <- round(seq(min(node.size.raw), max(node.size.raw), length.out = 5), digits = 1)
            axis(side = 4, at = seq(0, 1, length.out = 5), labels = n.legend.labels, pos = 1, cex.axis = legend.axis.cex)

            edge_legend <- as.raster(matrix(rev(raw.edge.colors), ncol=1))
            plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pair-Score',
                 cex.main = legend.title.cex)
            rasterImage(edge_legend, 0.75, 0, 1, 1)
            e.legend.labels <- round(seq(min(edge.widths), max(edge.widths), length.out = 5), digits = 1)
            axis(side = 4, at = seq(0, 1, length.out = 5), labels = e.legend.labels, pos = 1, cex.axis = legend.axis.cex)

        } else {

            plot(network, layout = coords, asp = 0, ...)

        }


    # otherwise, return igraph object
    } else {

        return(network)

    }

}
