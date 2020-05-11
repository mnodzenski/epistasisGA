#' A function to plot a network of potential epistatic snps
#'
#' This function rplot a network of potential epistatic snps.
#'
#' @param results.df The \code{unique.results} data frame of results of GA runs from \code{combine.islands}, with fitness scores restricted to the top results using function \code{network.threshold}.
#' @param chrom.size An integer indicating the number of snps in the chromosome.
#' @param seed an integer specifying the graph seed
#' @return A list whose first element is the fitness score and second element is the sum of weighted difference vectors for the target snps.
#'
#' @importFrom  dplyr group_by summarize %>%
#' @import igraph
#' @export

network.plot <- function(results.df, chrom.size, seed = 10){

  set.seed(seed)
  n.top.chroms <- nrow(results.df)
  result.df$h.score <- results.df$fitness.score*results.df$n.islands.found

  all.edge.weights <- do.call(rbind, lapply(1:n.top.chroms, function(res.row){

    chrom.res <- results.df[res.row, ]
    chrom <- as.vector(t(chrom.res[ , 1:chrom.size]))
    chrom.pairs <- expand.grid(chrom, chrom)
    chrom.pairs <- chrom.pairs[chrom.pairs$Var1 != chrom.pairs$Var2, ]
    chrom.pair.edges <- do.call(rbind, lapply(1:nrow(chrom.pairs), function(pair.row){

      pair <- chrom.pairs[pair.row, ]
      pair1 <- pair[ , 1]
      pair2 <- pair[ , 2]
      if (pair1 > pair2){

          pair[ , 1] <- pair2
          pair[ , 2] <- pair1
          pair1 <- pair[ , 1]
          pair2 <- pair[ , 2]

      }

      #pair1.pos <- which(chrom.res[ , 1:chrom.size] == pair1)
      #pair2.pos <- which(chrom.res[ , 1:chrom.size] == pair2)
      #standardized.elements <- as.vector(t(chrom.res[ , (1+chrom.size):(2*chrom.size)]))
      #vec.length.sq <- sum(standardized.elements^2)
      #pair1.val <- standardized.elements[pair1.pos]
      #pair2.val <- standardized.elements[pair2.pos]
      #pair.contrib <- sum(standardized.elements[c(pair1.pos, pair2.pos)]^2)
      #pair.weight <- (pair.contrib/vec.length.sq)*chrom.res$fitness.score
      pair.weight <- chrom.res$h.score
      pair$edge.weight <-  pair.weight
      #pair$var1.val <- pair1.val
      #pair$var2.val <- pair2.val
      return(pair)

    }))

    return(chrom.pair.edges)

  }))

  final.edge.weights <- all.edge.weights %>%
                          group_by(Var1, Var2) %>%
                          summarize(edge.weight = 1/sum(edge.weight)) %>%
                          as.data.frame()

  node.scores <- data.frame(name = unlist(results.df[1:n.top.chroms, 1:chrom.size]),
                            h.score = unlist(results.df[1:n.top.chroms, "h.score"]))

  sum.node.scores <- node.scores %>%
    group_by(name) %>%
    summarize(size = sum(h.score)) %>%
    as.data.frame()

  colnames(final.edge.weights)[1:2] <- c("from", "to")
  network <- graph.data.frame(final.edge.weights[ , 1:2], directed = F, vertices = sum.node.scores)
  E(network)$weight <- final.edge.weights$edge.weight
  color_fun <- colorRampPalette(c('red', 'grey', 'white'))
  required.colors <- as.integer(as.factor(E(network)$weight))
  colors <- color_fun(length(required.colors))
  E(network)$color <- colors[required.colors]

  color_fun <- colorRampPalette(c('white', 'grey', 'green'))
  required.colors <- as.integer(as.factor(V(network)$size))
  colors <- color_fun(length(required.colors))
  V(network)$color <- colors[required.colors]
  V(network)$shape <- "crectangle"

  plot(network, layout = layout.kamada.kawai, vertex.label.cex = 0.5, asp = 0)
}
