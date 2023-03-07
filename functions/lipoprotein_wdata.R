lipoprotein_wdata = function(wd = df, define_Dstate = "fasting"){
  wdata = wd %>% filter( grepl( "Lipoprotein subclasses", class) & dietary_state == define_Dstate ) %>%
  dplyr::select( outcome, subclass, raw.label, Obs_beta, Obs_P, MR_beta, MR_P)
  ##
  wdata$trait = "concentration"
  ##
  w = grep("Total lipids in", wdata$raw.label);
  q = grep("%", wdata$raw.label[w]); if(length(q)>0){w = w[-q]};  wdata$trait[w] = "total_lipids"
  ##
  w = grep("Phospholipids in ", wdata$raw.label)
  q = grep("%", wdata$raw.label[w]); if(length(q)>0){w = w[-q]}; wdata$trait[w] = "phospholipids"
  ##
  w = grep("Total cholesterol in", wdata$raw.label); 
  q = grep("%", wdata$raw.label[w]); if(length(q)>0){w = w[-q]}; wdata$trait[w] = "total_cholesterol"
  ##
  w = grep("Cholesterol esters in", wdata$raw.label); 
  q = grep("%", wdata$raw.label[w]); if(length(q)>0){w = w[-q]}; wdata$trait[w] = "cholesterol_esters"
  ##
  w = grep("Free cholesterol in", wdata$raw.label); 
  q = grep("%", wdata$raw.label[w]); if(length(q)>0){w = w[-q]}; wdata$trait[w] = "free_cholesterol"
  ##
  w = grep("Triglycerides in", wdata$raw.label); 
  q = grep("%", wdata$raw.label[w]); if(length(q)>0){w = w[-q]}; wdata$trait[w] = "triglycerides"
  #if(define_Dstate != "response"){
    ## RATIO
    w = grep("Phospholipids in ", wdata$raw.label)
    q = grep("%", wdata$raw.label[w]); if(length(q)>0){w = w[q]}; wdata$trait[w] = "phospholipids/total_lipids"
    ## RATIO
    w = grep("Total cholesterol in", wdata$raw.label); 
    q = grep("%", wdata$raw.label[w]); if(length(q)>0){w = w[q]}; wdata$trait[w] = "total_cholesterol/total_lipids"
    ## RATIO
    w = grep("Cholesterol esters in", wdata$raw.label); 
    q = grep("%", wdata$raw.label[w]); if(length(q)>0){w = w[q]}; wdata$trait[w] = "cholesterol_esters/total_lipids"
    ## RATIO
    w = grep("Cholesterol esters in", wdata$raw.label); 
    q = grep("%", wdata$raw.label[w]); if(length(q)>0){w = w[q]}; wdata$trait[w] = "cholesterol_esters/total_lipids"
    ## RATIO
    w = grep("Free cholesterol in", wdata$raw.label); 
    q = grep("%", wdata$raw.label[w]); if(length(q)>0){w = w[q]}; wdata$trait[w] = "free_cholesterol/total_lipids"
    ## RATIO
    w = grep("Triglycerides in", wdata$raw.label); 
    q = grep("%", wdata$raw.label[w]); if(length(q)>0){w = w[q]}; wdata$trait[w] = "triglycerides/total_lipids"
    
    ## A vector to describe [concentrations] from ratios (%)
    wdata$value = "concentration"
    w = grep("/", wdata$trait); if(length(w)>0){ wdata$value[w] = "ratio" }
    
    ## remove "raito" from subclass names
    wdata$subclass = gsub(" ratios", "" ,wdata$subclass)
    
    ## Define Trait Order
    wdata$trait = factor(wdata$trait, levels = c("concentration", "total_lipids", 
                                               "phospholipids", 
                                               "total_cholesterol", "cholesterol_esters", "free_cholesterol",
                                               "triglycerides",
                                               ###
                                               "phospholipids/total_lipids", 
                                               "total_cholesterol/total_lipids", 
                                               "cholesterol_esters/total_lipids", 
                                               "free_cholesterol/total_lipids",
                                               "triglycerides/total_lipids")[length(unique(wdata$trait)):1] )
  # } else {
  #   ## Define Trait Order
  #   wdata$trait = factor(wdata$trait, levels = c("concentration", "total_lipids", 
  #                                              "phospholipids", 
  #                                              "total_cholesterol", "cholesterol_esters", "free_cholesterol",
  #                                              "triglycerides")[length(unique(wdata$trait)):1] )
  #   
  # }
  
  ## Define plot order by defining factor levels
  wdata$subclass = factor(wdata$subclass, levels = c("Extremely large VLDL","Very large VLDL",
                                                                 "Large VLDL","Medium VLDL","Small VLDL","Very Small VLDL",
                                                                 "IDL",
                                                                 "Large LDL","Medium LDL","Small LDL",
                                                                 "Very large HDL", "Large HDL","Medium HDL","Small HDL"))
  
  
  ## return data to user
  return(wdata)
}
