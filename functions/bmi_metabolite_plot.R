bmi_metabolite_plot = function(data, 
                           metabolite, 
                           metabolite_name = NA,
                           color_choices = NA){
  ## define plot colors
  if( !is.na(color_choices) ){
    pcol = RColorBrewer::brewer.pal(8, "Set1")  
  } else {
    pcol = color_choices
  }
  
  ## define meatbolites
  fasted = paste0(metabolite, "_f")
  postprandial = paste0(metabolite, "_p")
  response = paste0(metabolite, "_r")
  
  ## define the plot data
  plot_data = na.omit( data[, c("sex", "age", "bmi", fasted, postprandial, response )] )
  
  ## remove outlires
  w = unique(c( id_outliers(plot_data[, fasted] ), 
         id_outliers(plot_data[, postprandial] ), 
         id_outliers(plot_data[, response] ) ) )
  if(length(w)>0){ plot_data = plot_data[-w,] }
  
  ## Define BMI Class
  plot_data = as.data.frame(plot_data)
  plot_data$bmi_class = "healthy weight"
  w = which( plot_data$bmi > 25  )
  plot_data$bmi_class[w] = "overweight"
  w = which( plot_data$bmi > 30 )
  plot_data$bmi_class[w] = "obese"
  w = which( plot_data$bmi > 40 )
  plot_data$bmi_class[w] = "severely obese"
  
  plot_data$bmi_class = factor(plot_data$bmi_class, 
                               levels = c("healthy weight","overweight","obese", "severely obese") )
  
  
  ## FASTING PLOT
  library(mgcv)
  form = formula( paste0(fasted, " ~ sex + s(age) + bmi") )
  fit0 = gam(form  , data = plot_data)
  form = formula( paste0(fasted, " ~ sex + s(age) + s(bmi)") )
  fit = gam(form, data = plot_data)
  lrt = lmtest::lrtest(fit0, fit)
  if(lrt$LogLik[1] == lrt$LogLik[2]){
    lrt_p = 1
  } else { 
    # lrt$Df[2] < 1 
    lrt_p = lrt$"Pr(>Chisq)"[2]
    }
  ##
  p1 = plot_data %>% ggplot(aes_string(x = "bmi", y = fasted) ) +
    geom_point( color = "grey80" , alpha = .50, size = 3) +
    geom_smooth( method = lm, formula = y~x, color = "blue", fill = "blue") +
    geom_smooth( method = gam, formula = y~s(x), 
                 color = "black", fill = "black",
                 size = 3, se = FALSE ) +
    geom_smooth( method = lm, formula = y~x, aes(color = bmi_class)  ) +
    scale_color_brewer(palette="Set1") +
    labs(color = "BMI category", x = "BMI", y = paste0("fasting ", metabolite) ) +
    ggtitle( paste0("gam lrt P = ", formatC(lrt_p) ) ) +
    theme_bw()
  
  ## POSTPRANDIAL PLOT
  form = formula( paste0(postprandial, " ~ sex + s(age) + bmi") )
  fit0 = gam(form  , data = plot_data)
  form = formula( paste0(postprandial, " ~ sex + s(age) + s(bmi)") )
  fit = gam(form, data = plot_data)
  lrt = lmtest::lrtest(fit0, fit)
  if( lrt$Df[2] < 1 ){
    lrt_p = 1
  } else { lrt_p = lrt$"Pr(>Chisq)"[2] }
  ##
  p2 = plot_data %>% ggplot(aes_string(x = "bmi", y = postprandial)) +
    geom_point( color = "grey80" , alpha = .50, size = 3) +
    geom_smooth( method = lm, formula = y~x, color = "blue", fill = "blue") +
    geom_smooth( method = gam, formula = y~s(x), 
                 color = "black", fill = "black",
                 size = 3, se = FALSE ) +
    geom_smooth( method = lm, formula = y~x, aes(color = bmi_class)  ) +
    scale_color_brewer(palette="Set1") +
    labs(color = "BMI category", x = "BMI", y = paste0("postprandial ", metabolite)  ) +
    ggtitle( paste0("gam lrt P = ", formatC(lrt_p) ) ) +
    theme_bw()
  
  ## FASTING V POSTPRANDIAL PLOT
  p3 = plot_data %>% ggplot(aes_string(x = fasted, y = postprandial)) +
    geom_point( color = "grey80" , alpha = .50, size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", size = 1.25) +
    #geom_vline(xintercept = c(0.3, 0.4, 0.5), linetype = "dotted") +
    #geom_hline(yintercept = c(0.3, 0.4, 0.5, 0.6), linetype = "dotted") +
    geom_smooth( method = lm, formula = y~x-1, color = "black" ) +
    geom_smooth( method = lm, formula = y~x-1, aes(color = bmi_class)  ) +
    scale_color_brewer(palette="Set1") +
    theme_bw() +
    labs(color = "BMI category", x = paste0("fasting ", metabolite), y = paste0("postprandial ", metabolite) )
  
  ## RESPONSE PLOT
  form = formula( paste0(response, " ~ sex + s(age) + bmi") )
  fit0 = gam(form  , data = plot_data)
  form = formula( paste0(response, " ~ sex + s(age) + s(bmi)") )
  fit = gam(form, data = plot_data)
  lrt = lmtest::lrtest(fit0, fit)
  if( lrt$Df[2] < 1 ){
    lrt_p = 1
  } else { lrt_p = lrt$"Pr(>Chisq)"[2] }
  ##
  p4 = plot_data %>% ggplot(aes_string(x = "bmi", y = response)) +
    geom_point( color = "grey80" , alpha = .50, size = 3) +
    geom_smooth( method = lm, formula = y~x, color = "blue", fill = "blue") +
    geom_smooth( method = gam, formula = y~s(x), 
                 color = "black", fill = "black",
                 size = 3, se = FALSE ) +
    geom_smooth( method = lm, formula = y~x, aes(color = bmi_class)  ) +
    scale_color_brewer(palette="Set1") +
    labs(color = "BMI category", x = "BMI", y = paste0(metabolite, " response") ) +
    ggtitle( paste0("gam lrt P = ", formatC(lrt_p) ) ) +
    theme_bw()
  
  
  #PLOT = ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, 
  #                 common.legend = TRUE, 
  #                labels = c("A", "B", "C", "D"))
  library(patchwork)
  PLOT = (p1|p2) / (p3|p4) & theme(legend.position = "bottom")
  PLOT = PLOT + 
    plot_annotation(title = metabolite_name, tag_levels = 'A') + 
    plot_layout(guides = "collect")
  
  return(PLOT)
  
}
