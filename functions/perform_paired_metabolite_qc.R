paired_metabolite_qc = function(wdata, 
                             trait1,
                             trait2,
                             single_trait_IQR_distance = 10, 
                             paired_delta_distance = 5){
  ## Define data frame
  df = wdata[, c(trait1, trait2)]
  rownames(df) = paste0("sample_", 1:nrow(df))
  
  ## turn NAs (from metaboprep) into zeros
  df[is.na(df)] = 0
  ##############################
  ## IDENTIFY the zeros
  ##############################
  ## Identify the zero values
  trait1_zero = which( df[,1] == 0 )
  trait2_zero = which( df[,2] == 0 )
  w = unique(c(trait1_zero, trait2_zero))
  df$zeros = "nonzero"
  if(length(w)>0){ df$zeros[w] = "zero" }
  ## define factor and levels
  df$zeros = factor(df$zeros, levels = c("nonzero","zero"))
  
  ##############################
  ## identify extreme outliers
  ##############################
  ## ID outliers
  trait1_iqr = id_outliers( df[-trait1_zero,1], 
                            outlier_method = "iqr", outlier_cutoff = single_trait_IQR_distance)
  trait2_iqr = id_outliers( df[-trait2_zero,2] , 
                            outlier_method = "iqr", outlier_cutoff = single_trait_IQR_distance)
  
  ## Add outlier ID to data frame 'df'
  df$iqr = as.character(df$zeros)
  df$iqr[ df$iqr == "nonzero" ] = "good"
  if(length(trait1_iqr)>0){ df$iqr[-trait1_zero][trait1_iqr] = "iqr_outlier" }
  if(length(trait2_iqr)>0){ df$iqr[-trait2_zero][trait2_iqr] = "iqr_outlier" }
  ## define factor and levels
  df$iqr = factor(df$iqr, levels = c("good","zero","iqr_outlier"))
  
  ##############################
  ## identify  outliers based on
  ## the two timepoints
  ##############################
  ## combined vector of positions
  zero_iqr_removal = which(df$iqr != "good")
  samples_2_analyze = which(df$iqr == "good")
  
  ###################
  ## (2) Delta
  ###################
  ## estimate a euclidean delta
  df$delta_outliers = as.character( df$iqr )
  
  if(length(zero_iqr_removal)>0){
    xydelta = apply(df[samples_2_analyze,1:2], 1, function(z){ z[2]-z[1] })
    xydelta_outliers = id_outliers( xydelta, 
                                 outlier_method = "iqr", outlier_cutoff = paired_delta_distance)
    if(length(xydelta_outliers)>0){ df$delta_outliers[samples_2_analyze][xydelta_outliers] = "delta_outlier" }
  } else {
    xydelta = apply(df[,1:2], 1, function(z){ z[2]-z[1] })
    xydelta_outliers = id_outliers( xydelta, 
                                 outlier_method = "iqr", outlier_cutoff = paired_delta_distance)
    if(length(xydelta_outliers)>0){ df$delta_outliers[xydelta_outliers] = "delta_outlier" }
  }
  ## define factor and levels
  df$delta_outliers = factor(df$delta_outliers, levels = c("good","zero","iqr_outlier", "delta_outlier"))
  
  ## Return data frame
  return(df)
}
