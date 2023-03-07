dh_forrest_plot = function(data, 
                           alpha = 0.05, 
                           analysis_type = "tsls",
                           arrange_by = "pval",
                           title = NA, 
                           covariates = NA, 
                           confounders = NA,
                           outcome_label = "",
                           exposure_label = "",
                           filter_data = TRUE,
                           column_4_x = "analysis",
                           column_4_color = "analysis",
                           col_count = 1,
                           color_choices = NA,
                           facet_variable = "outcome",
                           facet_bg_label_color = "grey",
                           facet_text_label_color = "black"){
  
  if(filter_data == TRUE){
    ## identify the data to plot by alpha (p-value) and analysis focus
    outcomes2plot = as.character( unlist( tibble(data) %>%
                                            filter(pval <= alpha & analysis %in% analysis_type ) %>%
                                            arrange_at(arrange_by) %>% dplyr::select(outcome) ) )
    
    ## define plot data
    ## define pdata; used from making plot
    w = which(data$outcome %in% outcomes2plot)
    plot_data = tibble( data[w, ] )
    
    ## Set metabolite plotting order
    plot_data$outcome = factor(plot_data$outcome, levels = outcomes2plot)
  }else{
    plot_data = tibble( data )
  }
  
  #################
  ## plot variables
  #################
  rowcount = nrow(plot_data)
  l = length( unique( unlist(plot_data[, column_4_color] )) )
  if( is.na(color_choices) ){
    pcol = RColorBrewer::brewer.pal(l, "Set1")[l:1]
  } else {
    pcol = color_choices[l:1]
  }
  
  
  #################
  ## plot
  #################
  PLOT = plot_data %>% ggplot( aes_string(x = column_4_x, y = "beta",
                                          ymin = "lowerCI" , ymax = "upperCI" ) ) +
    #geom_pointrange( aes( col = analysis , shape = sig ) , size = 1) +
    geom_pointrange( aes_string( col = column_4_color , shape = "sig" ) , size = 0.5) +
    # scale_size(range = c(0.5, 1.5)) +
    scale_colour_manual( values = pcol ) +
    scale_shape_manual( values = c( not_sig = 21, sig = 19), drop = FALSE ) +
    geom_hline( yintercept = 0, linetype = 2 ) +
    xlab(outcome_label) +
    ylab(exposure_label) +
    geom_errorbar( aes_string(ymin="lowerCI", ymax="upperCI", col=column_4_color),
                   width = 0.75 , cex = 0.5 ) +
    theme_bw() +
    facet_wrap(as.formula(paste("~", facet_variable)),  ncol = col_count, strip.position="left", nrow = rowcount, scales = "free_y") +
    theme( axis.text.y=element_blank(), axis.ticks.y=element_blank(),
           strip.text.y = element_text( size = 10, color = facet_text_label_color )) +
    theme(strip.background = element_rect(fill = facet_bg_label_color ) )+
    coord_flip() +
    labs(shape = "significance") +
    guides(colour = guide_legend(reverse = T)) +
    theme(strip.text.y.left = element_text(angle = 0))
  
  if( !is.na(title) ){
    PLOT = PLOT + labs( title  = title )
  }
  
  if( !is.na(title) & !is.na(covariates) ){
    PLOT = PLOT + labs( title  = title,
                        subtitle = paste0("covariates: ", paste0(covariates, collapse = ", ")  ) )
  }
  
  if( !is.na(title) & !is.na(covariates) & !is.na(confounders) ){
    PLOT = PLOT + labs( title  = title,
                        subtitle = paste0("covariates: ", paste0(covariates, collapse = ", ") ,
                                          "\nconfounders: " , paste0(confounders, collapse = ", ") ) )
  }
  
  return(PLOT)
  
}
