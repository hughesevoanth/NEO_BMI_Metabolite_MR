paired_outlier_filtering = function(wdata, 
                             trait1,
                             trait2,
                             IQR_distance = 5, 
                             paired_SD_distance = 10, 
                             paired_delta_distance = 10, 
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
  
  ##############################
  ## identify extreme outliers
  ##############################
  ## ID outliers
  trait1_iqr = id_outliers( df[-trait1_zero,1], 
                            outlier_method = "iqr", outlier_cutoff = IQR_distance)
  trait2_iqr = id_outliers( df[-trait2_zero,2] , 
                            outlier_method = "iqr", outlier_cutoff = IQR_distance)
  
  ## Add outlier ID to data frame 'df'
  df$iqr = "good"
  if(length(trait1_iqr)>0){ df$iqr[-trait1_zero][trait1_iqr] = "outlier" }
  if(length(trait2_iqr)>0){ df$iqr[-trait2_zero][trait2_iqr] = "outlier" }
  
  
  ##############################
  ## identify  outliers based on
  ## the two timepoints
  ##############################
  ## remove zeros and IQR outliers from the following QC steps
  a = which(df$zeros == "zero")
  b = which(df$iqr == "outlier")
  zero_iqr_removal = unique(c(a,b))
  df$zero_iqr_removal = "good"
  if(length(zero_iqr_removal)>0){
    df$zero_iqr_removal[zero_iqr_removal] = "outliers"
  }
  ###################
  ## (1) Standard Deviations
  ###################
  ## estimate the standard deviation between the two values
  df$sd_outliers = NA
  
  if(length(zero_iqr_removal)>0){
    xySD = apply(df[-zero_iqr_removal,1:2], 1, function(z){ sd(z, na.rm = TRUE) })
    SD = median(xySD, na.rm = TRUE)
    ## Define the outliers
    df$sd_outliers[-zero_iqr_removal] = apply(df[-zero_iqr_removal,1:2], 1, function(z){
            ifelse( abs(z[2]-z[1]) <= SD*paired_SD_distance, "good" , "outlier") 
          })
  } else {
    xySD = apply(df[,1:2], 1, function(z){ sd(z, na.rm = TRUE) })
    SD = median(xySD, na.rm = TRUE)
    ## Define the outliers
    df$sd_outliers = apply(df[,1:2], 1, function(z){
            ifelse( abs(z[2]-z[1]) <= SD*paired_SD_distance, "good" , "outlier") 
          })
  }
  
  ###################
  ## (2) Delta
  ###################
  ## estimate a euclidean delta
  df$delta_outliers = NA
  
  if(length(zero_iqr_removal)>0){
    xydelta = apply(df[-zero_iqr_removal,1:2], 1, function(z){ abs(z[2]-z[1]) })
    xydelta = apply(df[-zero_iqr_removal,1:2], 1, function(z){ z[2]-z[1] })
    DELTA = median(xydelta, na.rm = TRUE)
    cutoff = DELTA*paired_delta_distance
    
    ## Define the outliers
    df$delta_outliers[-zero_iqr_removal] = sapply(xydelta,  function(z){
          ifelse( z <= cutoff, "good" , "outlier") 
        })
  } else {
    xydelta = apply(df[,1:2], 1, function(z){ abs(z[2]-z[1]) })
    DELTA = median(xydelta, na.rm = TRUE)
    cutoff = DELTA*paired_delta_distance
    
    ## Define the outliers
    df$delta_outliers = sapply(xydelta,  function(z){
          ifelse( z <= cutoff, "good" , "outlier") 
        })
  }
  
  
  ###################
  ## (3) linear model residuals
  ###################
  ## build a linear model
  if(length(zero_iqr_removal)>0){ 
      temp = df[-zero_iqr_removal, 1:2] 
  } else {
      temp = df[, 1:2]
    }
  colnames(temp) = c("fast","post")
  fit = lm(post ~ fast, data = temp)
  # res = residuals(fit)
  fitQ = rq(post ~ fast, data = temp)
  res = residuals(fitQ)
  m = match(rownames(temp), names(res))
  ##
  temp$res = res[m]
  temp$cooksdist = cooks.distance(fit)[m]
  mean_res = mean(temp$res, na.rm = TRUE)
  sd_res = sd(temp$res, na.rm = TRUE)
  cutoff = mean_res + (sd_res * residual_distance)
  
  ## Define the residual outliers
  temp$residual_outliers = sapply(temp$res,  function(z){
          ifelse( abs(z) <= cutoff, "good" , "outlier") 
        })
  
  ## Define the cooks distance outliers
  temp$cooks_outliers = sapply(temp$cooksdist,  function(z){
          ifelse( z <= cooks_distance, "good" , "outlier") 
        })
  
  ## add linear model data to data frame df
  m = match(rownames(df), rownames(temp))
  df = cbind(df, temp[m,-c(1:2)])
  
  ## Return data frame
  return(df)
}
