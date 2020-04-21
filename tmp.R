pathway <-  examplePathways[["5991130_Programmed_Cell_Death"]]
stats <- exampleRanks
gseaParam <- 1

rnk <- rank(-stats)
ord <- order(rnk)
statsAdj <- stats[ord]
statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
statsAdj <- statsAdj/max(abs(statsAdj))
pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
pathway <- sort(pathway)
gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                        returnAllExtremes = TRUE)
bottoms <- gseaRes$bottoms
tops <- gseaRes$tops
n <- length(statsAdj)
xs <- as.vector(rbind(pathway - 1, pathway))
ys <- as.vector(rbind(bottoms, tops))
toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
diff <- (max(tops) - min(bottoms))/8
x = y = NULL
g <- ggplot(toPlot, aes(x = x, y = y)) +
  geom_point(color = "green",
             size = 0.1) + geom_hline(yintercept = max(tops), colour = "red",
                                      linetype = "dashed") +
  geom_hline(yintercept = min(bottoms),colour = "red", linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "black") +
  geom_line(color = "green") +
  theme_bw() +
  geom_segment(data = data.frame(x = pathway),
               mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), size = ticksSize) +
  theme(panel.border = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "rank", y = "enrichment score")


creat_plot_data <- function(query_idx, ordered_vec){

  query_idx <- sort(query_idx)
  es <-
    calcGseaStat(ordered_vec, query_idx, returnAllExtremes = TRUE)
  n <- length(ordered_vec)
  xs <- as.vector(rbind(query_idx - 1, query_idx))
  ys <- as.vector(rbind(es$bottoms, es$tops))
  toPlot <-
    data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))

  return(toPlot)

}

up_genes <- c("G3", "G4", "G8")
down_genes <- c("G12", "G13", "G14", "G19")

up_gene_idx <- match(up_genes, names(test_sig))
down_gene_idx <- match(up_genes, names(test_sig))

es_up <- calcGseaStat(test_sig, up_gene_idx,returnAllExtremes = TRUE, returnLeadingEdge = TRUE)
es_down <- calcGseaStat(test_sig, down_gene_idx)


