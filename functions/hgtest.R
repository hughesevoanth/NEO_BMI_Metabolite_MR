hgtest = function(allfeatures, selectedfeatures) {
  ## unique list of categories in the selected feature list
  cats = na.omit( unique(selectedfeatures) )
  ## iterate over those cat and perform test
  
  HG_F_test = sapply(cats, function(i){
    g = length(allfeatures)
    ## total number of tested features with some annotation
    N = length(allfeatures) - ( sum(is.na(allfeatures) | allfeatures == "") )
    ## number of features in cat i
    m = sum(allfeatures == i, na.rm = TRUE )
    ## number of features NOT in cat i
    n = N - m
    ## number of features selected by X criteria, that have some annotation
    k = length(selectedfeatures) - ( sum(is.na(selectedfeatures) | selectedfeatures == "") )
    ## number of selected features in cat i
    x = sum(selectedfeatures == i, na.rm = TRUE)
    
    
    ## estimate fold enrichment and p.value
    fold.enrichment <- (x / k ) / (m / N)
    p.value <- phyper(q=x-1, m=m, n=n, k=k, lower.tail=FALSE)
    
    ## FISHER EXACT TEST
    ## BUILD 2 X 2 contingency table
    dmat = matrix( c(x,k-x,m-x,n-(k-x)),2,2, byrow=TRUE, dimnames = list(c("Selected","NotSelected"), c("FocusedCat", "AllOtherCats") ))
    #print(i)
    #print(dmat)
    ftest <- fisher.test(x=dmat, alternative = "greater")
    
    
    ## data out
    out = c(fold.enrichment, p.value, ftest$estimate, ftest$p.value)
    names(out) = c("fold.enrichment", "HG_pval", "Ftest_OR", "Ftest_pval")
    return(out)
  })
  ## return to user
  return(t(HG_F_test))
  
  
}

