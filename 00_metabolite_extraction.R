######################################
## Metabolite Extraction
##  - NEO BMI Metabolite MR -
##  Prepare data to run in metaboprep
## by: David Hughes
## date: Aug 18 2021
######################################

#########################
## READ IN THE DATA
#########################
### parameter file
pfile = "parameters/pfile.txt"
par_data = read.table(pfile, header = FALSE, sep = "=", as.is = TRUE)

### Read in the metabolite data
w = grep("pheno_file", par_data[,1])
f = paste0(par_data[w,2] )
pdata = read.table(f, header = TRUE, sep = ",", as.is = TRUE, fill = FALSE)
pdata$ID = paste0("ind", pdata$ID)

####################
## Remove the "." in 
## column names
####################
## the code here is not just a simple gsub(".","_",x), as it 
## was creating some errors

k = sapply(colnames(pdata), function(i){
  substring(i, nchar(i)-1, nchar(i))
} )
w = which(k == ".1")
if(length(w)>0){
  colnames(pdata)[w] = sapply( colnames(pdata)[w], function(x){
    paste0( substring(x, 1, nchar(x)-2 ) , "_1")
  })
}
#######
#######
k = sapply(colnames(pdata), function(i){
  substring(i, nchar(i)-1, nchar(i))
} )
w = which(k == ".3")
if(length(w)>0){
  colnames(pdata)[w] = sapply( colnames(pdata)[w], function(x){
    paste0( substring(x, 1, nchar(x)-2 ) , "_3")
  })
}


#############################
##
## DEFINE METABOLITE TRAITS
##
#############################

### FASTING DATA
fast = grep("_1", colnames(pdata))[-c(1:2)]
fast = colnames(pdata)[fast]

### POSTPRANDIAL DATA
post = grep("_3", colnames(pdata))[-c(1:2)]
post = colnames(pdata)[post]

### RESPONSE
respo = grep( "r", substr(colnames(pdata), 1,1) )[-1]
respo = colnames(pdata)[respo]

## Trait Vector
traits = c(fast, post, respo)


#############################
##
## Extract metabolite data
##
#############################
metdata = pdata[, traits]
rownames(metdata) = paste0("ind", rownames(metdata))

batchdata = pdata[, c("neolight","visitdd")]
rownames(batchdata) = rownames(metdata)

#############################
##
## Annotate metabolite data
##
#############################
## annotation file from the metaboprep package
ng_anno = metaboprep:::ng_anno
feature_data = data.frame(ids = traits)
####################
## edit trait names
####################
id2 = traits
## remove time stamp 1 and 3
id2 = gsub( "_1","", id2)
id2 = gsub( "_3","", id2)
## remove timestamp "r"
id2 = sapply(id2, function(x){
  first = substring(x, 1,1)
  if(first == "r"){
    out = substring(x, 2,nchar(x))
  } else {
    out = x
  }
  return(out)
})
## define as character
id2 = as.character(id2)
## change the percentages "p" to "pct"
id2 = sapply(id2, function(x){
  last = substring(x, nchar(x), nchar(x))
  if(last == "p" & nchar(x)>3 ){
    out = paste0(x, "ct")
  } else {
    out = x
  }
  return(out)
})
## define as character
id2 = as.character(id2)
## make lowercase
id2 = tolower(id2)
## add to data frame
feature_data$metabolite =  id2 

####################
## define time point
####################
feature_data$dietary_state = "response"
w = grep("_1",feature_data$ids)
feature_data$dietary_state[w] = "fasting"
w = grep("_3",feature_data$ids)
feature_data$dietary_state[w] = "postprandial"

####################
## make a new ID for each metabolite
####################
newid = sapply(1:nrow(feature_data), function(i){
  state = feature_data$dietary_state[i]
  met = feature_data$metabolite[i]
  if(state == "fasting"){
    met = paste0(met, "_f")
  }
  if(state == "postprandial"){
    met = paste0(met, "_p")
  }
  if(state == "response"){
    met = paste0(met, "_r")
  }
  return(met)
})

feature_data$new_id = newid
  
####################
## match names
####################
m = match(feature_data$metabolite, ng_anno$metabolite)
feature_data = cbind(feature_data, ng_anno[m,-1] )

#############################
##
## Write to file
##
#############################
w = grep("data_dir", par_data[,1] )
filepath = paste0(par_data[w,2], "source/release2/")
## make a metaboprep data folder in data source directory
cmd = paste0("mkdir -p ", filepath, "metaboprep")
system(cmd)
### write metabolite data to file
f =  paste0(filepath, "metaboprep/source_metabolites_extracted.txt")
write.table(metdata, file = f, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

### write batchdata data to file
f =  paste0(filepath, "metaboprep/source_metabolite_batch_data.txt")
write.table(batchdata, file = f, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

### write feature data to file
f =  paste0(filepath, "metaboprep/source_metabolites_feature_data.txt")
write.table(feature_data, file = f, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)



