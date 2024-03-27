################################
## Build NEO BMI-GRS
## by: David Hughes
## date: January 7th 2021
################################
library(tidyverse)
library(readxl)

pfile = "parameters/pfile.txt"
pdata = read.table(pfile, header = FALSE, sep = "=", as.is = TRUE)

## Read in the genotype Data
## geno files
genofiles = c("source/locke_SNPs.csv","source/Yengo_SNPs.csv")
genodata = lapply(genofiles, function(file){
	w = grep("data_dir", pdata[,1])
	f = paste0(pdata[w,2], file)
	out = read.table(f, header = TRUE, sep = ",", as.is = TRUE)
	})
dosage = list( locke = t(genodata[[1]]) , yengo = t(genodata[[2]])  )
colnames(dosage[[1]]) = dosage[[1]][1,]
colnames(dosage[[2]]) = dosage[[2]][1,]
dosage[[1]] = dosage[[1]][-1,]
dosage[[2]] = dosage[[2]][-1,]
##
rownames(dosage[[1]]) = gsub( "exm.rs","rs",rownames(dosage[[1]]) )
rownames(dosage[[2]]) = gsub( "exm.rs","rs",rownames(dosage[[2]]) )


## Read in the GWAS files with the beta effect estimates
gwasfiles = c("Locke_BMI_SNPs_TableS9_HetAncestry.txt","yengo_BMI_656.txt")
gwasdata = lapply(gwasfiles, function(file){
	w = grep("data_dir", pdata[,1])
	f = paste0(pdata[w,2], file)
	out = read.table(f, header = TRUE, sep = "\t", as.is = TRUE)
	})
colnames(gwasdata[[1]])[8] = "BETA"
## rs9581854 is also known as rs12016871
w = which(gwasdata[[1]]$SNP == "rs12016871")
gwasdata[[1]]$SNP[w] = "rs9581854"


## Read in the SNP info file
snpfiles = c("locke_SNP_info.xlsx","Yengo_SNP_info.xlsx")
snpinfo = lapply(snpfiles, function(file){
	w = grep("data_dir", pdata[,1])
	f = paste0(pdata[w,2],"source/" ,file)
	out = read_excel(f)
	})

### edit rsids LOCKE
snpinfo[[1]]$rsid = gsub("exm.","",snpinfo[[1]]$SNP)
w = which( !is.na( snpinfo[[1]]$Note) )
snpinfo[[1]]$rsid[w] = snpinfo[[1]]$Note[w]
### edit rsids YENGO
snpinfo[[2]]$rsid = snpinfo[[2]]$'Note (SNP names in extraction file, exchanged names used in Exome chip)'
snpinfo[[2]]$rsid = gsub("exm.","",snpinfo[[2]]$rsid)
### edit SNP id
snpinfo[[1]]$SNP = gsub("-",".",snpinfo[[1]]$SNP)
snpinfo[[2]]$SNP = gsub("-",".",snpinfo[[2]]$SNP)


######################################
##
## Add SNP INFO to dosage data 
##
######################################
## LOCKE
m = match( rownames(dosage[[1]]), snpinfo[[1]]$rsid )
dosage[[1]] = cbind( snpinfo[[1]][m,], dosage[[1]] )

## YENGO
m = match( rownames(dosage[[2]]), snpinfo[[2]]$rsid )
dosage[[2]] = cbind( snpinfo[[2]][m,], dosage[[2]] )

######################################
##  **** LOCKE *****
## Add GWAS effect estimate data to 
## dosage file
##
######################################
## rs9581854 is also known as rs12016871
######################################
new_data = c()
####################
for(i in 1:nrow(dosage[[1]]) ){
	## define SNP ID
	snpid = dosage[[1]]$rsid[i]
	## define effect allele
	ea = dosage[[1]]$'Effect allele'[i]
	aa = dosage[[1]]$'alternative allele'[i]
	## define observed EAF
	eaf = dosage[[1]]$'EAF'[i]
	########################
	## find SNP in gwas data
	########################
	w = which( gwasdata[[1]]$SNP == snpid )
	if(length(w) == 1){
		## define gwas SNP data
		gwas_ea = gwasdata[[1]]$effect_allele[w]
		gwas_aa = gwasdata[[1]]$alt_allele[w]
		## define gwas EAF
		gwas_eaf = gwasdata[[1]]$EAF[w]
		## Do the alleles match ?
		allele_match = sum(ea == gwas_ea & aa == gwas_aa)
		## do the EAF match?
		eaf_match = abs(eaf - gwas_eaf)
		if( allele_match == 1 & eaf_match <= 0.1 ){
			## LOCKE
			new_data = rbind(new_data, gwasdata[[1]][w, c("Chr","Position","BETA") ])
			## YENGO
			# new_data = rbind(new_data, gwasdata[[1]][w, c("CHR","POS","BETA_COJO") ])
		} else {
			allele_flip = sum(ea == gwas_aa & aa == gwas_ea)
			eaf_flip = abs(eaf - (1-gwas_eaf) )
			if(allele_flip == 1 & eaf_flip <= 0.1 ){
				out = gwasdata[[1]][w, c("Chr","Position","BETA") ]
				out[3] = out[3] * -1
				new_data = rbind( new_data, out )
			} 
		}
	} else {
				new_data = rbind( new_data, c("Chr","Position","NA") )
			}
}
####### redefine rownames
rownames(new_data) = rownames(dosage[[1]])
#######  re-define dosage
dosage[[1]] = cbind(dosage[[1]]$SNP, new_data, dosage[[1]][,-1])
colnames(dosage[[1]])[1] = "SNP"


####################################
## if any effect estimates are 
##  negative then flip the alleles
####################################
locke = dosage[[1]]
w = which( locke$BETA < 0)
##
dos = which(colnames(locke) %in% 1:5744)
##
for(i in w){
	locke$BETA[i] = locke$BETA[i]  * -1
	###
	locke[i, dos] = 2 - as.numeric( locke[i, dos] )
	##
	ea = locke$'Effect allele'[i]
	aa = locke$'alternative allele'[i]
	eaf = locke$'EAF'[i]
	##
	locke$'Effect allele'[i] = aa
	locke$'alternative allele'[i] = ea
	locke$'EAF'[i] = 1-as.numeric(eaf)
}


## write LOCKE data to file
w = grep("data_dir", pdata[,1])
f = paste0(pdata[w,2], "processed/locke_dosage_data_QCd.txt")
write.table(locke, file = f, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

# locke = read.table(f, header = TRUE, sep = "\t", as.is = TRUE)

######################################
##  **** YENGO *****
## Add GWAS effect estimate data to 
## dosage file
##
######################################
## rs9581854 is also known as rs12016871
######################################
new_data = c()
####################
for(i in 1:nrow(dosage[[2]]) ){
	## define SNP ID
	snpid = dosage[[2]]$rsid[i]
	## define effect allele
	ea = dosage[[2]]$'EA'[i]
	aa = dosage[[2]]$'NEA'[i]
	## define observed EAF
	eaf = dosage[[2]]$'EAF'[i]
	########################
	## find SNP in gwas data
	########################
	w = which( gwasdata[[2]]$SNP == snpid )
	if(length(w) == 1){
		## define gwas SNP data
		gwas_ea = gwasdata[[2]]$Tested_Allele[w]
		gwas_aa = gwasdata[[2]]$Other_Allele[w]
		## define gwas EAF
		gwas_eaf = gwasdata[[2]]$Freq_Tested_Allele_in_HRS[w]
		## Do the alleles match ?
		allele_match = sum(ea == gwas_ea & aa == gwas_aa)
		## do the EAF match?
		eaf_match = abs(eaf - gwas_eaf)
		if( allele_match == 1 & eaf_match <= 0.1 ){
			## LOCKE
			# new_data = rbind(new_data, gwasdata[[1]][w, c("Chr","Position","BETA") ])
			## YENGO
			new_data = rbind(new_data, c(snpid, gwasdata[[2]][w, c("CHR","POS","BETA") ] ) )
		} else {
			allele_flip = sum(ea == gwas_aa & aa == gwas_ea)
			eaf_flip = abs(eaf - (1-gwas_eaf) )
			if(allele_flip == 1 & eaf_flip <= 0.1 ){
				out = c(snpid, gwasdata[[2]][w, c("CHR","POS","BETA") ] )
				out[4] = as.numeric(out[4]) * -1
				new_data = rbind( new_data, out )
			} else {
				new_data = rbind( new_data, c(snpid, "Chr","Position","NA") )
			} 
		}
	} else {
				new_data = rbind( new_data, c(snpid, "Chr","Position","NA") )
			}
}
####### redefine rownames
rownames(new_data) = new_data[,1]
new_data = new_data[, -1]
new_data = as.data.frame(new_data)

for(i in 1:ncol(new_data)){
	new_data[,i] = unlist(new_data[,i])
}

####
dosage[[2]] = cbind( dosage[[2]]$SNP, new_data, dosage[[2]][,-1])
colnames(dosage[[2]])[1] = "SNP"
# colnames(dosage[[2]])[4] = "BETA"

####################################
## if any effect estimates are 
##  negative then flip the alleles
####################################
yengo = dosage[[2]]
w = which(is.na( as.numeric( unlist( yengo$BETA)) ) )
yengo = yengo[-w, ]
##
w = which( as.numeric( unlist( yengo$BETA ) ) < 0)
##
dos = which(colnames(yengo) %in% 1:5744)
##
for(i in w){
	yengo$BETA[i] = as.numeric( yengo$BETA[i] ) * -1
	###
	yengo[i, dos] = 2 - as.numeric( yengo[i, dos] )
	##
	ea = yengo$'EA'[i]
	aa = yengo$'NEA'[i]
	eaf = yengo$'EAF'[i]
	##
	yengo$'EA'[i] = aa
	yengo$'NEA'[i] = ea
	yengo$'EAF'[i] = 1-as.numeric(eaf)
}

####
colnames(yengo)[10] = "note"

#yengo$CHR = unlist(yengo$CHR)
#yengo$POS = unlist(yengo$POS)
#yengo$BETA = unlist(yengo$BETA)

## write LOCKE data to file
w = grep("data_dir", pdata[,1])
f = paste0(pdata[w,2], "processed/yengo_dosage_data_QCd.txt")
write.table(yengo, file = f, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

# yengo = read.table(f, header = TRUE, sep = "\t", as.is = TRUE)




####################################
## 
## Generate GRS
##
####################################
#############
## LOCKE
#############
dos = which(colnames(locke) %in% 1:5744)
locke_grs = t(apply(locke[, dos], 2, function(x){
	un_w_grs = sum( as.numeric(x) , na.rm = TRUE)
	w_grs = sum( as.numeric(x) * as.numeric(locke$BETA) , na.rm = TRUE)
	out = c(un_w_grs, w_grs)
	return(out) 
	}))
###
colnames(locke_grs) = c("Locke_Un_Weighted_GRS","Locke_Weighted_GRS")
#############
## YENGO
#############
dos = which(colnames(yengo) %in% 1:5744)
yengo_grs = t(apply(yengo[, dos], 2, function(x){
	un_w_grs = sum( as.numeric(x) , na.rm = TRUE)
	w_grs = sum( as.numeric(x) * as.numeric(yengo$BETA) , na.rm = TRUE)
	out = c(un_w_grs, w_grs)
	return(out) 
	}))
###
colnames(yengo_grs) = c("Yengo_Un_Weighted_GRS","Yengo_Weighted_GRS")

###
grs = cbind(locke_grs, yengo_grs)

## WRITE TO FILE
w = grep("data_dir", pdata[,1])
f = paste0(pdata[w,2], "processed/NEO_BMI_GRS.txt")
write.table(grs, file = f, row.names = TRUE,
	col.names = TRUE, sep = "\t", quote = FALSE)

## correlation between the two GRS's
cor.test(grs[,1] , grs[,2]) ## r = 0.914
cor.test(grs[,1] , grs[,3]) ## 0.3520288
cor.test(grs[,2] , grs[,4]) ## 0.5737177
cor.test(grs[,3] , grs[,4]) ## 0.8939391





