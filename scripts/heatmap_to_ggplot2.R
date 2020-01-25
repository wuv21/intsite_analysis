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

heatmap_to_ggplot2 <- function(heatmap_rds, heatmap_type, metric_order, output_suffix, sep_levels) {
  stopifnot(heatmap_type %in% c("genomic", "epigenetic"))
  
  rds <- readRDS(heatmap_rds)
  
  roc <- rds[["ROC"]]
  pval <- rds$pvalues$np
  
  if (heatmap_type == "epigenetic") {
    rownames(roc) <- gsub("\\.10Kb", "", rownames(roc))
    rownames(pval) <- gsub("\\.10Kb", "", rownames(pval))
    
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
  
  pval <- data.frame(pval[rownames(pval) %in% metric_order, ],
                     check.names = FALSE)
  
  pval$y <- rownames(pval)
  pval <- pval %>%
    gather(sampleName, np_val, -y) %>%
    mutate(y = factor(y, levels = rev(metric_order)),
           np_val = as.numeric(np_val))
  
  df <- data.frame(roc[rownames(roc) %in% metric_order, ],
                   check.names = FALSE)
  
  df$y <- rownames(df)
  df <- df %>%
    gather(sampleName, roc, -y) %>%
    mutate(y = factor(y, levels = rev(metric_order)))
  
  df <- left_join(df, pval, by = c("y", "sampleName")) %>%
    mutate(pval_lbl = ifelse(np_val < 0.001, "***",
                             ifelse(np_val < 0.01, "**",
                                    ifelse(np_val < 0.05, "*", ""))))
  
  graph_settings <- list(
    geom_tile(aes(fill = roc), color = "#222222", size = 0.5),
    geom_text(aes(label = pval_lbl), size = 7),
    scale_fill_gradient2(low = scale_low,
                         high = scale_hi,
                         mid = "white",
                         midpoint = 0.5,
                         limits = c(0,1)),
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 0.5),
          axis.text.y = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text =  element_text(size = 15),
          legend.position = "none",
          panel.background = element_blank(),
          panel.grid = element_blank()),
    coord_equal())
  
  if (heatmap_type == "genomic") {
    g <- ggplot(df, aes(x = sampleName, y = y)) +
      graph_settings
    
    p1_x <- ggplot_build(g)$layout$panel_params[[1]]$x.labels
    p1_y <- ggplot_build(g)$layout$panel_params[[1]]$y.labels
    
    p <- egg::set_panel_size(g,
                             height=unit(length(p1_y), "cm"),
                             width=unit(length(unique(p1_x)), "cm"))
    
    p2 <- make_xaxis_plot(p1_x, sep_levels, xaxis_colors, sep = "_")
    
    p_comb <- cowplot::plot_grid(p, p2, ncol = 1, align = "v")
    
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
    
    p1_x <- ggplot_build(g1)$layout$panel_params[[1]]$x.labels
    p2_x <- ggplot_build(g2)$layout$panel_params[[1]]$x.labels
    
    s_height_1 <- unit(length(unique(df1$y)), "cm")
    s_width_1 <- unit(length(unique(df1$sampleName)), "cm")
    
    s_height_2 <- unit(length(unique(df2$y)), "cm")
    s_width_2 <- unit(length(unique(df2$sampleName)), "cm")
    
    g1 <- egg::set_panel_size(g1, height = s_height_1, width = s_width_1)
    g2 <- egg::set_panel_size(g2, height = s_height_2, width = s_width_2)
    
    p1 <- make_xaxis_plot(p1_x, sep_levels, xaxis_colors, sep = "_")
    p1_g1 <- cowplot::plot_grid(g1, p1, ncol = 1, align = "v")
    
    p2 <- make_xaxis_plot(p2_x, sep_levels, xaxis_colors, sep = "_")
    p2_g2 <- cowplot::plot_grid(g2, p2, ncol = 1, align = "v")
    
  }
  
  general_fn <- glue("figs/{heatmap_type}_heatmap_{output_suffix}")
  save_height <- 24
  
  if (heatmap_type == "genomic") {
    cowplot::save_plot(filename = glue("{general_fn}.svg"),
                       plot = p_comb,
                       device = "svg",
                       base_height = save_height)
  } else {
    cowplot::save_plot(filename = glue("{general_fn}_1.svg"),
                       plot = p1_g1,
                       device = "svg",
                       base_height = save_height)
    
    cowplot::save_plot(filename = glue("{general_fn}_2.svg"),
                       plot = p2_g2,
                       device = "svg",
                       base_height = save_height)
  }
  
  return(TRUE)
}