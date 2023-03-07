lipid_conc_plot_by_bmi_cat = function(wdata = NULL, 
                           trait_vector = NULL, 
                           ylabel = NULL, 
                           scales = "fixed"){
  ########################
  ## Define BMI categories
  ########################
  cats = c("healthy","overweight","obese","severely obese")
  
  ########################
  ## Define the data
  ########################
  lipid_conc_ss = c()
  ##
  for(lipo in trait_vector){
    w = which( colnames(wdata) %in% lipo )
    d = wdata[,w]
    out = sapply(cats, function(cat){
      q = which(wdata$bmi_cat == cat)
      o = c( lipo, cat, mean(d[q], na.rm = TRUE), quantile(d[q], probs = c(0.25, 0.75), na.rm = TRUE) )
      names(o) = c("metabolite","bmi","mean","lci","uci")
      return(o)
    })
    out = t(out)
    out = as.data.frame(out)
    for(i in 3:5){out[,i] = as.numeric(out[,i])}
    lipid_conc_ss = rbind(lipid_conc_ss, out)
  }
  ## Define factor levels
  lipid_conc_ss$metabolite = factor(lipid_conc_ss$metabolite, levels = trait_vector)
  lipid_conc_ss$bmi = factor(lipid_conc_ss$bmi, levels = cats)
  
  
  ## Make the plot
  PLOT = lipid_conc_ss %>% ggplot(aes(x = bmi, y = mean)) +
    geom_point(aes(color = bmi), size = 3) +
    geom_errorbar(aes(ymin = lci, ymax = uci, color = bmi), linewidth = 0.25) +
    facet_wrap(.~metabolite, nrow = 1, scales = scales) +
    theme_bw() +
    theme( axis.text.x = element_blank() ) +
    theme(legend.position="bottom") +
    labs(x = "BMI category", y = ylabel , color = "BMI category")
  
  ## Return the plot
  return(PLOT)
  
}
