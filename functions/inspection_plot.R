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
                                paired_SD_distance = 5, 
                                paired_delta_distance = 5, 
                                residual_distance = 5,
                                cooks_distance = 0.025)
  ################################
  ## Set up the colors
  ################################
  set1cols = RColorBrewer::brewer.pal(8, "Set1")
  pcol = c("grey",set1cols[2],set1cols[3],set1cols[1])
  
  ################################
  ## Build the individual plots
  ################################
  plot1 = df %>% ggplot(aes_string(x = fasting, y = postprandial)) +
    geom_point(aes(color = sd_outliers), size = 1) +
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
    geom_point(aes(color = delta_outliers), size = 1) +
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
    geom_point(aes(color = res_outliers), size = 1) +
    scale_color_manual(values = pcol, drop = FALSE) +
    ## equivalency line
    geom_abline( intercept = 0, slope = 1, 
                 color = "grey30", 
                 linetype = "dashed"  ) +
    # ## linear
    # geom_smooth(data = df %>% filter(zero_iqr_removal == "good"),
    #             method = "lm", formula = y ~ x, color = "black", size = 0.5 ) +
    ## median
    geom_quantile(data = df %>% filter(iqr == "good"), 
                  quantiles = 0.5, color = "blue", size = 1,
                  formula = y ~ x) +
    # ## gam
    # geom_smooth(data = df %>% filter(zero_iqr_removal == "good"),
    #             method = "gam", formula = y ~ s(x), color = "green", size = 0.5 ) +
    labs(title = "Median Regression Residual Outliers", color = "QC step") +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(size = 2.5)))
  
  plot4 = df %>% ggplot(aes_string(x = fasting, y = postprandial)) +
    geom_point(aes(color = cooks_outliers), size = 1) +
    scale_color_manual(values = pcol, drop = FALSE) +
    ## equivalency line
    geom_abline( intercept = 0, slope = 1, 
                 color = "grey30", 
                 linetype = "dashed"  ) +
    ## linear
    geom_smooth(data = df %>% filter(iqr == "good"),
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
