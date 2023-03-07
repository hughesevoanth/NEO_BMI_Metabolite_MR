paired_outlier_filtering = function(wdata, 
                             trait1,
                             trait2,
                             IQR_distance = 10, 
                             paired_SD_distance = 5, 
                             paired_delta_distance = 5, 
                             residual_distance = 5,
                             cooks_distance = 0.025){
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
                            outlier_method = "iqr", outlier_cutoff = IQR_distance)
  trait2_iqr = id_outliers( df[-trait2_zero,2] , 
                            outlier_method = "iqr", outlier_cutoff = IQR_distance)
  
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
  ## (1) Standard Deviations
  ###################
  ## estimate the standard deviation between the two values
  # df$sd_outliers = NA
  df$sd_outliers = as.character( df$iqr )
  
  if(length(zero_iqr_removal)>0){
    xySD = apply(df[samples_2_analyze,1:2], 1, function(z){ sd(z, na.rm = TRUE) })
    xySD_outliers = id_outliers( xySD, 
                 outlier_method = "iqr", outlier_cutoff = paired_SD_distance)
    if(length(xySD_outliers)>0){ df$sd_outliers[samples_2_analyze][xySD_outliers] = "sd_outlier" }
    } else {
      xySD = apply(df[,1:2], 1, function(z){ sd(z, na.rm = TRUE) })
      xySD_outliers = id_outliers( xySD, 
                                   outlier_method = "iqr", outlier_cutoff = paired_SD_distance)
      if(length(xySD_outliers)>0){ df$sd_outliers[xySD_outliers] = "sd_outlier" }
  }
  ## define factor and levels
  df$sd_outliers = factor(df$sd_outliers, levels = c("good","zero","iqr_outlier", "sd_outlier"))
  
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
  
  ###################
  ## (3) median regression 
  ##     and 
  ##     linear model residuals
  ###################
  temp = df[samples_2_analyze, 1:2] 
  colnames(temp) = c("fast","post")
  
  ## Median Regression
  fitQ = rq(post ~ fast, data = temp)
  ## Residuals
  res = residuals(fitQ)
  ## Match and add to temp data frame
  m = match(rownames(temp), names(res))
  temp$res = res[m]
  ## Identify median residual outliers
  res_outliers = id_outliers( temp$res, 
                                  outlier_method = "iqr", outlier_cutoff = residual_distance)
  temp$res_outliers = "good"
  if(length(res_outliers)>0){ temp$res_outliers[res_outliers] = "res_outlier" }
  
  ## Define linear model cooks distance outliers
  fit = lm(post ~ fast, data = temp)
  temp$cooksdist = cooks.distance(fit)[m]
  temp$cooks_outliers = sapply(temp$cooksdist,  function(z){
          ifelse( z <= cooks_distance, "good" , "cook_outlier") 
        })
  
  ## add linear model data to data frame df
  m = match(rownames(df), rownames(temp))
  df = cbind(df, temp[m,-c(1:2)])
  
  ## add zero and iqr outliers to res_outlier and cooks_outlier
  w = which(as.character(df$iqr) != "good")
  if(length(w)>0){
    df$res_outliers[w] = as.character( df$iqr[w] )
    df$cooks_outliers[w] = as.character( df$iqr[w] )
  }
  ## define factor and levels
  df$res_outliers = factor(df$res_outliers, levels = c("good","zero","iqr_outlier", "res_outlier"))
  df$cooks_outliers = factor(df$cooks_outliers, levels = c("good","zero","iqr_outlier", "cook_outlier"))
  
  ## Return data frame
  return(df)
}
