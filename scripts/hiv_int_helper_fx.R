# make hg38 reference file
make_ref <- function(file) {
  # make reference file
  if (!file.exists(file)) {
    ref <- load_ref_files(list(file = "refSeqComposite:refGene", symbol = "name2"),
                          type = "GRanges",
                          freeze = "hg38")
    
    saveRDS(ref, file = "reference_hg38.rds")
    
  } else {
    ref <- readRDS("reference_hg38.rds")
    
    # change ref listData
    og_names <- names(ref@elementMetadata@listData)
    og_names[13] <- "geneName"
    names(ref@elementMetadata@listData) <- og_names
    
    char_col <- c("geneName", "name", "chrom", "strand", "geneName")
    for (i in char_col) {
      ref@elementMetadata@listData[[i]] <- as.character(ref@elementMetadata@listData[[i]])
    }
    
    rm_col <- c("bin", "name", "exonStarts", "exonEnds", "score", "cdsStartStat", "cdsEndStat", "exonFrames", "annot_sym")
    for (i in rm_col) {
      ref@elementMetadata@listData[[i]] <- NULL
    }
  }
  
  return(ref)
}

# calc_abundance
calc_abund <- function(primer, posid) {
  if (length(posid) == 1) {
    return(1)
  }
  
  ids <- list()
  for (i in seq_along(posid)) {
    current_idp <- posid[i]
    
    if (length(ids) == 0) {
      ids[[current_idp]] <- list()
      
      if (primer[i] == "u3") {
        ids[[current_idp]][[1]] <- 1
        ids[[current_idp]][[2]] <- 0
        
      } else {
        ids[[current_idp]][[1]] <- 0
        ids[[current_idp]][[2]] <- 1       
      }
      
      next
    }
    
    if (current_idp %in% names(ids)) {
      if (primer[i] == "u3") {
        ids[[current_idp]][[1]] <- ids[[current_idp]][[1]] + 1
        
      } else {
        ids[[current_idp]][[2]] <- ids[[current_idp]][[2]] + 1   
      }
      
    } else {
      ids[[current_idp]] <- list()
      
      if (primer[i] == "u3") {
        ids[[current_idp]][[1]] <- 1
        ids[[current_idp]][[2]] <- 0
        
      } else {
        ids[[current_idp]][[1]] <- 0
        ids[[current_idp]][[2]] <- 1       
      }
    }
  }
  
  if (length(ids) > 2) {
    print("WARNING")
    print(posid)
  }
  
  # compile average
  abund_sum <- 0
  for (pos in ids) {
    if (pos[[1]] == 0) {
      abund_sum <- abund_sum + pos[[2]]
    } else if (pos[[2]] == 0) {
      abund_sum <- abund_sum + pos[[1]]
    } else {
      abund_sum <- abund_sum + base::mean(c(pos[[1]], pos[[2]]))
    }
    
  }
  return(abund_sum)
}

# find breakpoint
find_bp <- function(i, df) {
  row <- df[i, ]
  
  row$primer <- tolower(row$primer)
  
  # find middle of breakpoint (+2bp from edge)
  pos <- -1
  ori <- ""
  if (row$primer == "u5" & row$strand == "+") {
    pos <- row$start + 2
    ori <- "+"
  } else if (row$primer == "u5" & row$strand == "-") {
    pos <- row$end - 2
    ori <- "-"
  } else if (row$primer == "u3" & row$strand == "+") {
    pos <- row$start + 2
    ori <- "-"
  } else {
    pos <- row$end - 2
    ori <- "+"
  }
  
  pos_full <- paste(row$seqnames, ori, pos, sep = "")
  
  return(pos_full)
}

# draw chord diagram
draw_chord <- function(df_chord, uniq_cells, cell_colors, abund_scale) {
  circos.clear()
  factors = 1:length(df_chord$posID)
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  
  circos.par(cell.padding = c(0, 0, 0, 0),
             gap.degree = 360/(length(df_chord$posID) + 1),
             canvas.xlim = c(-1.75, 1.75), canvas.xlim = c(-1.75, 1.75))
  
  circos.initialize(factors, xlim = c(0, 25))
  
  # for sample_cell
  circos.track(ylim = c(0, 1),
               track.height = 0.1,
               bg.border = df_chord$cell_track_color)
  
  # sites
  circos.track(ylim = c(1, 2),
               track.height = 0.1,
               bg.border = df_chord$abund_color)
  
  uniq <- !duplicated(df_chord$sample_cell)  ## logical vector of unique values 
  uniq_idx <- seq_along(df_chord$sample_cell)[uniq]  ## indices 
  uniq_val <- df_chord$sample_cell[uniq] ## the values 
  
  for (i in c(1: (length(uniq_idx)))) {
    p1 <- uniq_idx[i]
    p2 <- -1
    
    if (i == length(uniq_idx)) {
      p2 <- length(df_chord$sample_cell)
    } else if (p2 == p1 + 1) {
      p2 <- p1
    } else {
      p2 <- uniq_idx[i + 1] - 1
    }
    
    # highlight.sector(c(p1: p2),
    #                  track.index = 1,
    #                  col = paste(cell_colors[which(uniq_val[i] == uniq_cells)], "60", sep = ""),
    #                  text = uniq_val[i],
    #                  facing = "downward",
    #                  niceFacing = TRUE,
    #                  text.vjust = "120mm",
    #                  text.col = cell_colors[which(uniq_val[i] == uniq_cells)],
    #                  cex = 3)

    highlight.sector(c(p1: p2),
                     track.index = 1,
                     col = paste(cell_colors[which(uniq_val[i] == uniq_cells)], "60", sep = ""))
  }
  
  df_intersects <- df_chord %>%
    dplyr::filter(value != 0)
  
  for (i in seq_along(df_intersects$posID)) {
    df_row <- df_intersects[i, ]
    
    se <- c(df_row$from_posID, df_row$uniq_posID)
    circos.link(se[1],
                0,
                se[2],
                0,
                lwd = 2,
                col = "#000000")
    
    
    highlight.sector(c(df_row$from_posID: df_row$from_posID),
                     track.index = 1,
                     text = df_row$gene_plot,
                     col = "#FFFFFFFF",
                     facing = "clockwise",
                     # niceFacing = TRUE,
                     text.vjust = "45mm",
                     cex = 1.5)
  }
  
  # label more abundant clones
  # df_abundant <- df_chord %>%
  #   dplyr::filter(abund >= 2)
  # 
  # for (i in seq_along(df_abundant$gene_plot)) {
  #   df_row <- df_abundant[i, ]
  #   
  #   highlight.sector(c(df_row$from_posID: df_row$from_posID),
  #                    track.index = 1,
  #                    text = df_row$gene_plot_noid,
  #                    col = "#FFFFFFFF",
  #                    facing = "clockwise",
  #                    niceFacing = TRUE,
  #                    text.vjust = "40mm",
  #                    cex = 1.5)
  # }
  
  circos.clear()
}

# rbind dfs
df_rbind <- function(df, v, cols) {
  df_tmp <- data.frame(v)
  colnames(df_tmp) <- cols
  
  if (is.null(dim(df))) {
    df <- df_tmp
  } else {
    df <- rbind(df, df_tmp)
  }
  
  return(df)
}

# pairwise stat pooling
pairwise_stat_pool <- function(df_pair, p) {
  subset_comb <- combn(unique(df_pair$cell_sample), 2)
  
  df_sim <- NULL
  for(i in c(1:dim(subset_comb)[2])) {
    df_pair_sub <- df_pair %>%
      filter(cell_sample %in% subset_comb[, i])
    
    ids <- df_pair_sub$posID
    weights <- df_pair_sub$abund
    
    actual_intersect <- duplicated(df_pair_sub$posID) | duplicated(df_pair_sub[nrow(df_pair_sub):1, "posID"])[nrow(df_pair_sub):1]
    actual_intersect_sum <- length(unique(df_pair_sub[actual_intersect, "posID"]))
    
    
    pool_ids <- rep(ids, times = floor(weights))
    
    sub_1_bins <- df_pair_sub %>%
      filter(cell_sample == subset_comb[1, i]) %>%
      select(abund) %>%
      mutate(abund = floor(abund)) %>%
      sum(.)
    
    sub_2_bins <- df_pair_sub %>%
      filter(cell_sample == subset_comb[2, i]) %>%
      select(abund) %>%
      mutate(abund = floor(abund)) %>%
      sum(.)
    
    test <- c()
    for (j in c(1:10000)) {
      if (sub_1_bins + sub_2_bins != length(pool_ids)) {
        print(subset_comb[, i])
      }
      
      pool_scrambled <- sample(pool_ids)
      t1 <- pool_scrambled[c(1:sub_1_bins)]
      t2 <- pool_scrambled[c(sub_1_bins + 1 : length(pool_scrambled))]
      
      test <- c(test, length(intersect(t1, t2)))
    }
    
    tmp <- data.frame(overlaps = test,
                      sub1 = subset_comb[1, i],
                      sub2 = subset_comb[2, i],
                      act_sum = actual_intersect_sum)
    
    df_sim <- df_rbind(df_sim, tmp, colnames(tmp))
  }
  
  df_sim <- df_sim %>%
    mutate(sub_combo = paste(sub1, sub2, sep = " | "),
           sub_combo = factor(sub_combo))
  
  df_sim_act <- df_sim %>%
    group_by(sub_combo) %>%
    summarise(act = max(act_sum))
  
  g <- ggplot(df_sim, aes(x = overlaps)) +
    geom_density(color = "#cccccc", fill = "#dddddd", alpha = 0.5, adjust = 3) +
    theme_bw() +
    labs(title = "Random overlap simulations",
         subtitle = paste("10000 trials for patient ", p, " | blue line indicates actual overlap", sep = ""),
         x = "Number of overlapping integration sites",
         y = "Density") +
    expand_limits(x = 0, y = 0) +
    facet_wrap(~ sub_combo, scales = "free")
    
  g <- g +
    geom_vline(data = df_sim_act, aes(xintercept = act), color = "#0008ff")
  
  return(g)
}