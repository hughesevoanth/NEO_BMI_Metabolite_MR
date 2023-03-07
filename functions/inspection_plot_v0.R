inspection_plot = function(wdata, metabolite){
  ################################
  ## Define the to metabolite states
  ################################
  fasting = paste0(metabolite, "_f")
  postprandial = paste0(metabolite, "_p")
  ################################
  ## Identify the outliers
  ################################
  df = paired_outlier_filtering(wdata = wdata, 
                                trait1 = fasting,
                                trait2 = postprandial,
                                IQR_distance = 10, 
                                paired_SD_distance = 10, 
                                paired_delta_distance = 10, 
                                residual_distance = 5,
                                cooks_distance = 0.025)
  ################################
  ## Set up the colors
  ################################
  df$color = "good"
  w = which(df$zeros == "zero")
  df$color[w] = "zero"
  w = which(df$iqr == "outlier")
  df$color[w] = "iqr"
  
  ### SD outliers
  df$color_sd = df$color
  w = which(df$sd_outliers == "outlier")
  df$color_sd[w] = "sd"
  df$color_sd = factor(df$color_sd, levels = c("good","zero","iqr", "sd"))
  ### Delta outliers
  df$color_delta = df$color
  w = which(df$delta_outliers == "outlier")
  df$color_delta[w] = "delta"
  df$color_delta = factor(df$color_delta, levels = c("good","zero","iqr", "delta"))
  ### Residual outlier
  df$color_residual = df$color
  w = which(df$residual_outliers == "outlier")
  df$color_residual[w] = "residual"
  df$color_residual = factor(df$color_residual, levels = c("good","zero","iqr", "residual"))
  ### Cook's distance outlier
  df$color_cooks = df$color
  w = which(df$cooks_outliers == "outlier")
  df$color_cooks[w] = "cooks"
  df$color_cooks = factor(df$color_cooks, levels = c("good","zero","iqr", "cooks"))
  
  set1cols = RColorBrewer::brewer.pal(8, "Set1")
  pcol = c("grey",set1cols[2],set1cols[3],set1cols[1])
  
  ################################
  ## Build the individual plots
  ################################
  plot1 = df %>% ggplot(aes_string(x = fasting, y = postprandial)) +
    geom_point(aes(color = color_sd), size = 1) +
    scale_color_manual(values = pcol, drop = FALSE) +
    ## equivalency line
    geom_abline( intercept = 0, slope = 1, 
                 color = "grey30", 
                 linetype = "dashed"  ) +
    # ## linear
    # geom_smooth(data = df %>% filter(zero_iqr_removal == "good"),
    #             method = "lm", formula = y ~ x, color = "black", size = 0.5 ) +
    # ## median
    # geom_quantile(data = df %>% filter(zero_iqr_removal == "good"), formula = y ~ x,
    #               quantiles = 0.5, color = "blue", size = 1) +
    # ## gam
    # geom_smooth(data = df %>% filter(zero_iqr_removal == "good"),
    #             method = "gam", formula = y ~ s(x), color = "green", size = 0.5 ) +
    ## title
    labs(title = "Standard Deviation Outliers", color = "QC step") +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(size = 2.5)))
  
  plot2 = df %>% ggplot(aes_string(x = fasting, y = postprandial)) +
    geom_point(aes(color = color_delta), size = 1) +
    scale_color_manual(values = pcol, drop = FALSE) +
    ## equivalency line
    geom_abline( intercept = 0, slope = 1, 
                 color = "grey30", 
                 linetype = "dashed"  ) +
    # ## linear
    # geom_smooth(data = df %>% filter(zero_iqr_removal == "good"),
    #             method = "lm", formula = y ~ x, color = "black", size = 0.5 ) +
    # ## median
    # geom_quantile(data = df %>% filter(zero_iqr_removal == "good"), 
    #               quantiles = 0.5, color = "blue", size = 1) +
    # ## gam
    # geom_smooth(data = df %>% filter(zero_iqr_removal == "good"),
    #             method = "gam", formula = y ~ s(x), color = "green", size = 0.5 ) +
    labs(title = "Fasting-Postprandial Delta Outliers", color = "QC step") +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(size = 2.5)))
  
  plot3 = df %>% ggplot(aes_string(x = fasting, y = postprandial)) +
    geom_point(aes(color = color_residual), size = 1) +
    scale_color_manual(values = pcol, drop = FALSE) +
    ## equivalency line
    geom_abline( intercept = 0, slope = 1, 
                 color = "grey30", 
                 linetype = "dashed"  ) +
    # ## linear
    # geom_smooth(data = df %>% filter(zero_iqr_removal == "good"),
    #             method = "lm", formula = y ~ x, color = "black", size = 0.5 ) +
    ## median
    geom_quantile(data = df %>% filter(zero_iqr_removal == "good"), 
                  quantiles = 0.5, color = "blue", size = 1) +
    # ## gam
    # geom_smooth(data = df %>% filter(zero_iqr_removal == "good"),
    #             method = "gam", formula = y ~ s(x), color = "green", size = 0.5 ) +
    labs(title = "Median Regression Residual Outliers", color = "QC step") +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(size = 2.5)))
  
  plot4 = df %>% ggplot(aes_string(x = fasting, y = postprandial)) +
    geom_point(aes(color = color_cooks), size = 1) +
    scale_color_manual(values = pcol, drop = FALSE) +
    ## equivalency line
    geom_abline( intercept = 0, slope = 1, 
                 color = "grey30", 
                 linetype = "dashed"  ) +
    ## linear
    geom_smooth(data = df %>% filter(zero_iqr_removal == "good"),
                method = "lm", formula = y ~ x, color = "black", size = 0.5 ) +
    # ## median
    # geom_quantile(data = df %>% filter(zero_iqr_removal == "good"), 
    #               quantiles = 0.5, color = "blue", size = 1) +
    # ## gam
    # geom_smooth(data = df %>% filter(zero_iqr_removal == "good"),
    #             method = "gam", formula = y ~ s(x), color = "green", size = 0.5 ) +
    labs(title = "Linear Model Cook's Distances Outliers", color = "QC step") +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(size = 2.5)))
  
  ################################
  ## patchwork plot
  ################################
  mytitle = paste0("Metabolite: ", metabolite)
  myplot = (plot1|plot2) / (plot3|plot4) +
    plot_annotation(title = mytitle, tag_levels = 'A' )
    # plot_layout(guides = "collect")
  
  #myplot = ggpubr::ggarrange(plot1, plot2, plot3, plot4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
  
  return(myplot)
}