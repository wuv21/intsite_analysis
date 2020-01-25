make_xaxis_plot <- function(g1_x, sep_levels, colors, sep, height_scale = 0.5) {
  df2 <- data.frame(org_x = factor(g1_x, levels = g1_x)) %>%
    separate(org_x, into = sep_levels, remove = FALSE, sep = sep) %>%
    pivot_longer(-org_x, names_to = "y", values_to = "y_val") %>%
    mutate(y = factor(y, levels = sep_levels))
  
  # g2_colors <- colors[unique(df2$y_val)]
  # 
  # print(names(g2_colors))
  # print(g2_colors)
  
  g2 <- ggplot(df2, aes(x = org_x, fill = y_val, y = reorder(y, plyr::desc(y)))) +
    geom_tile(color = "#000000") +
    theme_classic() +
    coord_equal() +
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line = element_blank()) +
    scale_fill_manual(values = colors)
  
  g2 <- egg::set_panel_size(g2,
                            height=unit(length(unique(df2$y)) * height_scale , "cm"),
                            width=unit(length(unique(g1_x)), "cm"))
  
  return(g2)
}