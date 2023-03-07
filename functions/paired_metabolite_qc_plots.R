paired_metabolite_qc_plots = function(wdata, metabolite){
  ################################
  ## Define the to metabolite states
  ################################
  fasting = paste0(metabolite, "_f")
  postprandial = paste0(metabolite, "_p")
  ################################
  ## Identify the outliers
  ################################
  df = paired_metabolite_qc(wdata = wdata, 
                            trait1 = fasting,
                            trait2 = postprandial,
                            single_trait_IQR_distance = 10, 
                            paired_delta_distance = 5)
  ################################
  ## Set up the colors
  ################################
  set1cols = RColorBrewer::brewer.pal(8, "Set1")
  pcol = c("grey",set1cols[2],set1cols[3],set1cols[1])
  
  ################################
  ## Deming Regression
  ################################
  form0 = as.formula( paste0(postprandial , " ~ ", fasting) )
  fitD = deming::deming(form0, data = df)
  
  ################################
  ## Build the individual plots
  ################################
  plot1 = df %>% ggplot(aes_string(x = fasting, y = postprandial)) +
    geom_point(aes(color = delta_outliers), size = 1) +
    scale_color_manual(values = pcol, drop = FALSE) +
    ## equivalency line
    geom_abline( intercept = 0, slope = 1, 
                 color = "grey30", 
                 linetype = "dashed"  ) +
    labs(title = "Raw data w/ outliers", color = "QC step") +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(size = 2.5)))
  
  plot2 = df %>% 
    filter(delta_outliers == "good") %>% 
    ggplot(aes_string(x = fasting, y = postprandial)) +
    geom_point(size = 1, color = "grey") +
    ## equivalency line
    geom_abline( aes(intercept = 0, slope = 1, color = "equivalency"), 
                 linetype = "dashed") +
    ## deming
    geom_abline( aes(intercept = fitD$coefficients[1], slope = fitD$coefficients[2], 
                     color = "Deming"),
                 size = 1.25) +
    # ## linear
    geom_smooth(method = "lm", formula = y ~ x, aes(color = "linear"), 
                size = 0.75, se = FALSE ) +
    # ## quadratic
    geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw=TRUE), 
                aes(color = "quadratic"), 
                size = 0.75, se = FALSE ) +
    # ## cubic
    geom_smooth(method = "lm", formula = y ~ poly(x, 3, raw=TRUE), 
                aes(color = "cubic"), 
                size = 0.75, se = FALSE ) +
    ## median
    geom_quantile(quantiles = 0.5, aes(color = "median"), size = 0.75,
                  formula = y ~ x) +
    ## gam
    geom_smooth(method = "gam", formula = y ~ s(x),
                aes(color = "gam"), se = FALSE, size = 0.75 ) +
    ## plot edits
    labs(title = "QC'd data") +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(size = 1.5))) +
    scale_color_manual(name="model", 
                       values = c( "equivalency" = "grey30",
                                   "Deming" = "black", 
                                   "linear" = "orangered1", 
                                   "quadratic" = "blue", 
                                   "cubic" = "purple", 
                                   "median" = "red3",
                                   "gam" = "green"),
                       guide = "legend") 
  
  
  ################################
  ## patchwork plot
  ################################
  mytitle = paste0("Metabolite: ", metabolite)
  myplot = plot1|plot2 + patchwork::plot_annotation(title = mytitle, tag_levels = 'A')
  
  return(myplot)
}
