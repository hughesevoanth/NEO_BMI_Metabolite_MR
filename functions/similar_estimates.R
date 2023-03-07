similar_estimates = function(x,y, scaledata = TRUE, sd_delta = 1){
  if(scaledata == TRUE){
    ## define data frame
    df = data.frame(x = scale(x), y = scale(y) )
    ##
    out1 = apply(df, 1, function(z){
      ifelse( abs(z[2]-z[1]) <= sd_delta , "same", "different")
      })
    out2 = apply(df, 1, function(z){
       ifelse( abs(z[1])>=sd_delta & abs(z[2])>=sd_delta, "extreme", "notextreme")
      })
    w = which(out1 == "same" & out2 == "extreme")
    if(length(w)>1){
      out1[w] = out2[w]
    }
    out = out1
    return(out)
  } else {
      ## define data frame
      df = data.frame(x = x, y = y )
      ##
      xySD = apply(df, 1, function(z){ sd(z, na.rm = TRUE) })
      SD = mean(xySD)
      SD = SD * sd_delta 
      out1 = apply(df, 1, function(z){
        ifelse( abs(z[2]-z[1]) <= SD, "same", "different"  )
        })
      out2 = apply(df, 1, function(z){
         ifelse( abs(z[1])>=SD & abs(z[2])>=SD, "extreme", "notextreme")
        })
      w = which(out1 == "same" & out2 == "extreme")
      if(length(w)>1){
        out1[w] = out2[w]
      }
      out = out1
      return(out)
    }
}