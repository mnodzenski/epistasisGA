#' A function to plot a network of potential epistatic snps
#'
#' This function rplot a network of potential epistatic snps.
#'
#' @param results.df The \code{unique.results} data frame of results of GA runs from \code{combine.islands}, with fitness scores restricted to the top results using function \code{network.threshold}.
#' @param repulse.rad A scalar affecting the graph shape. Decrease to reduce overlapping nodes.
#' @param node.size A scalar affecting the size of the graph nodes. Increse to increase size.
#' @param graph.area A scalar affecting the size of the graph area. Increase to increase graph area.
#' @param vertex.label.cex A scalar controlling the size of the vertex label. Increase to increase size.
#' @param seed an integer specifying the graph seed
#' @param plot A logical indicating whether the network should be plotted. If set to false, this function will return an igraph object which can be used for manual plotting.
#' @return An igraph object, if \code{plot} is set to FALSE.
#'
#' @importFrom  dplyr group_by summarize %>%
#' @import igraph
#' @importFrom qgraph qgraph.layout.fruchtermanreingold
#' @export

network.plot <- function(results.df, repulse.rad = 1000, node.size = 25, graph.area = 100,
                         vertex.label.cex = 0.5, seed = 10, plot = TRUE){

  set.seed(seed)
  chrom.size <- sum(grepl("snp", colnames(results.df)))/2
  n.top.chroms <- nrow(results.df)
  results.df$h.score <- results.df$fitness.score*results.df$n.islands.found
  results.df$h.score <- (results.df$h.score)/sd(results.df$h.score)

  all.edge.weights <- do.call(rbind, lapply(1:n.top.chroms, function(res.row){

    chrom.res <- results.df[res.row, ]
    chrom <- as.vector(t(chrom.res[ , 1:chrom.size]))
    chrom.pairs <- expand.grid(chrom, chrom)
    chrom.pairs <- chrom.pairs[chrom.pairs$Var1 != chrom.pairs$Var2, ]
    original.pair1 <- chrom.pairs$Var1
    original.pair2 <- chrom.pairs$Var2
    switch.these <- chrom.pairs$Var1 > chrom.pairs$Var2
    chrom.pairs$Var1[switch.these] <- original.pair2[switch.these]
    chrom.pairs$Var2[switch.these] <- original.pair1[switch.these]
    chrom.pairs <- chrom.pairs[!duplicated(chrom.pairs), ]
    chrom.pairs$h.score <- chrom.res$h.score
    chrom.pairs$fitness.score <- chrom.res$fitness.score
    return(chrom.pairs)

  }))

  final.edge.weights <- all.edge.weights %>%
                          group_by(Var1, Var2) %>%
                          summarize(edge.weight = 1/max(fitness.score)) %>%
                          as.data.frame()

  node.scores <- data.frame(name = unlist(results.df[1:n.top.chroms, 1:chrom.size]),
                            fitness.score = unlist(results.df[1:n.top.chroms, "fitness.score"]))

  node.scores <- node.scores %>%
    group_by(name) %>%
    summarize(size = max(fitness.score)) %>%
    as.data.frame()
  node.scores$size <- node.size*(node.scores$size/max(node.scores$size))
  st.weights <- (1/final.edge.weights$edge.weight)/max(1/final.edge.weights$edge.weight)

  colnames(final.edge.weights)[1:2] <- c("from", "to")
  network <- graph.data.frame(final.edge.weights[ , 1:2], directed = F, vertices = node.scores)
  E(network)$weight <- final.edge.weights$edge.weight
  E(network)$width <- st.weights
  color_fun <- colorRampPalette(c('red', 'grey', 'white'))
  required.colors <- as.integer(as.factor(E(network)$weight))
  colors <- color_fun(length(required.colors))
  edge.colors <- sapply(1:length(st.weights), function(x) adjustcolor(colors[required.colors][x], alpha.f = st.weights[x]))
  E(network)$color <- edge.colors

  color_fun <- colorRampPalette(c('white', 'grey', 'green'))
  required.colors <- as.integer(as.factor(V(network)$size))
  colors <- color_fun(length(required.colors))
  V(network)$color <- colors[required.colors]
  V(network)$shape <- "crectangle"

  if (plot){

    net.edges <- get.edgelist(network, names = F)
    coords <- qgraph.layout.fruchtermanreingold(net.edges, vcount = vcount(network),
                                                repulse.rad = repulse.rad*vcount(g),
                                                area = graph.area*(vcount(network)^2))
    plot(network, layout = coords, vertex.label.cex = vertex.label.cex, asp = 0)

  } else{

    return(network)

  }

}
