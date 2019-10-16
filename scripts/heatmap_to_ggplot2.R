pull_names <- function(query, prefix) {
  q_rownames <- rownames(query)
  matches <- startsWith(q_rownames, prefix)
  
  to_num <- gsub("k$", "000", q_rownames[matches])
  to_num <- gsub("M$", "000000", to_num)
  
  to_num <- str_sort(to_num, numeric = TRUE)
  
  to_num <- gsub("000000$", "M", to_num)
  to_num <- gsub("000$", "k", to_num)
  
  return(to_num)
}





heatmap_to_ggplot2 <- function(heatmap_rds, heatmap_type, metric_order, output_suffix) {
  stopifnot(heatmap_type %in% c("genomic", "epigenetic"))
  
  rds <- readRDS(heatmap_rds)
  
  roc <- rds[["ROC"]]
  
  if (heatmap_type == "epigenetic") {
    rownames(roc) <- gsub("\\.10Kb", "", rownames(roc))
    
    scale_low <- "purple"
    scale_hi <- "green"
  } else {
    metric_order <- c(pull_names(roc, "DNase"),
                      pull_names(roc, "CpG"),
                      pull_names(roc, "GC"),
                      pull_names(roc, "refSeq"),
                      "onco.100k",
                      "within_refSeq_gene",
                      "gene.width",
                      "general.width",
                      "start.dist",
                      "boundary.dist")
    
    scale_low <- "blue"
    scale_hi <- "red"
  }
  
  df <- data.frame(roc[rownames(roc) %in% metric_order, ],
                   check.names = FALSE)
  
  df$y <- rownames(df)
  df <- df %>%
    gather(sampleName, roc, -y) %>%
    mutate(y = factor(y, levels = rev(metric_order)))
  
  graph_settings <- list(
    geom_tile(aes(fill = roc), color = "#222222", size = 0.5),
    scale_fill_gradient2(low = scale_low,
                         high = scale_hi,
                         mid = "white",
                         midpoint = 0.5,
                         limits = c(0,1)),
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank()),
    coord_equal())
  
  if (heatmap_type == "genomic") {
    g <- ggplot(df, aes(x = sampleName, y = y)) +
      graph_settings
  } else {
    level_mid <- ceiling(length(levels(df$y)) / 2)
    
    df1 <- df %>%
      filter(as.numeric(y) > level_mid)
    
    df2 <- df %>%
      filter(as.numeric(y) <= level_mid)
    
    g1 <- ggplot(df1, aes(x = sampleName, y = y)) +
      graph_settings +
      theme(legend.position = "none")
    
    g2 <- ggplot(df2, aes(x = sampleName, y = y)) +
      graph_settings
    
    g <- grid.arrange(g1, g2, nrow = 1)
  }
  
  heatmap_png_fn <- paste0("figs/", heatmap_type, "_", "heatmap", "_", output_suffix, ".png")
  ggsave(filename = heatmap_png_fn, plot = g, device = "png", height = 10, width = 10)
  
  return(g)
}