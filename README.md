# NEO Metabolite MR

## -- fasting, postprandial, and response traits

	by: David Hughes
	date: Jan-May 2021
	

___

## Introduction

This repository containing the scripts used to carry out a one-sample MR analysis of the effects BMI on fasting, postprandial and response metabolite (Nightingale platform) traits in the The Netherlands Epidemiology of Obesity study (NEO).

## Metabolomics 2021

Please find in this repository a poster and video poster presentation for the Metabolomics 2021 conference.
				
## Scripts

### 01_generate_bmi_grs.R

This R sript was designed to be run inside the git repository directory where it can acces my hidden paramater folder and the paramater file contained within. 

This script:

	1. reads in the genotype dosage data for a set of requested SNPs from NEO
	2. reads in a file containing effect estimates (or weights) for each SNP on BMI
		i. Yengo et al estimates
		ii. Locke et al estimates 
	3. insures that the effect alleles are between dosage and effect file are aligned
	4. insures that all of the effect estimates are positive in orientation, otherwise the alleles are flipped
	5. then generates a GRS
		i. yengo un-weighted
		ii. yengo weighted
		iii. locke un-weighted
		iv. locke weighted
---

### 02_observational.Rmd

This is a Rmarkdown file that runs through observational analysis for BMI, metabolites, and covariates. The steps/analysis taken are:

	1. converting sampling date into a year month variable
	2. merging the phenotype and grs data together
	3. define a vector of confounders
		i. remove those with less than 20 observations
	4. evaluate the distribution of BMI and BMI-GRS
	5. compare distribution of BMI in NEO's high BMI sample and the Leiderdrop sample
	6. evaluate the sampling structure between NEO's high BMI and Leiderdrop sampling.
		i. do the genetic PCs correlate with the sampling populations?
			- yes :-( 
				- 3.3% of the variation in PC3 (p = 2.1x10-43)
				- 0.5% of the variation in PC4 (p = 4.35x10-8)
	7. evaluate distribution of shapiro W-statsitics across metbolites|dietary state.
	8. observational correlation between BMI and BMI-GRSs


|  GRS   					  | 	r 		 | 	 r2       |     pval     |
| ---------------------- | ---------- | --------- | ------------ |
|	Locke_Un_Weighted_GRS | 0.1105601 | 0.01222353 | 4.344514e-17 |
|	Locke_Weighted_GRS    | 0.1287896 | 0.01658675 | 1.130697e-22 |
|	Yengo_Un_Weighted_GRS | 0.2064412 | 0.04261798 | 2.521843e-56 |
|	Yengo_Weighted_GRS    | 0.2080231 | 0.04327361 | 3.501560e-57 |

	9. repeat observational correlations but with Rank normal transformed BMI.

|  GRS   					  | 	r 		 | 	 r2       |     pval     |
| ---------------------- | ---------- | --------- | ------------ |
|	Locke_Un_Weighted_GRS | 0.1187724 | 0.01410687 | 1.690487e-19 |
|	Locke_Weighted_GRS    | 0.1376140 | 0.01893761 | 1.098779e-25 |
|	Yengo_Un_Weighted_GRS | 0.2140008 | 0.04579635 | 1.737514e-60 |
|	Yengo_Weighted_GRS    | 0.2157269 | 0.04653811 | 1.848579e-61 |

	10. plot scatter plots of BMI and BMI-GRS relationships
	11. evaluate correlation between visit data (Year-Month) and BMI
		i. BMI (19.8% of variation) but not BMI-GRS (1.9%) is strongly influenced.
	12. evaluate effect of visit data (Year-Month) and BMI on metabolite abundance.
		i. univariate visit date on metabolite
		ii. univariate BMI on metabolite
		iii. multivariate visit date and BMI and metabolite
		iv. univariate BMI-GRS on metabolite
		v. multivariate visit date and BMI-GRS on metabolite
	13. Which metabolites are among the top influenced by visit date?
		- 88% of top 100 are response traits
	14. plot (boxplot) top two metabolites associated with visit date.
	15. How many of the top 100 metabolites associated with BMI are also among the top 100 metabolites associated with visit date? 
		- 0
	16. run iPVs (using spearman's rho and 1-abs(R) as the measure of distance) on metabolite data

|  tree cut height  | 	Me (number of PVs) | 
| ------------------| ----- | 
|	0.5 | 39 |
|	0.6 | 28 | 
|	0.7 | 16 | 

	17. estimate d.hat
	18. evaluate correlation between confounders and BMI, BMI-GRS, and dhat. 
		- whilst controling for sex and age
	19. evaluate correaltion between confounders and metabolites. 

---
### 03_tsls_mr.Rmd

This is a Rmarkdown file that runs through the one-sample MR. The steps/analysis taken are:


	1. converting sampling date into a year month variable
	2. merging the phenotype and grs data together
	3. define traits to run MR on
	4. run the observational and one sample MR [my function: onesample_mr() ]
	5. identify and table top observations [my functioni: top_onesample_mr() 
	6. convert onsample_mr() ouptput into long format [my function: convert_2_longformat_onesample_mr() ]
	7. generate forrest plot [my function: onesample_mr_forrest_plot() ]
	8. evaluate coefficient of variation across dietary states

---

## Description of functions

### onesample_mr()

	name: onesample_mr
	
	description: A function to perform observational and one-sample MR analysis.
		In all instances the outcome is rank normal transformed ( rntransform() ). 
		The function performs the following steps.
		1. observational analysis
			1. univariate analysis of exposure on outcome.
			2. multivariate analysis of exposure on outcome, with the inclusion of covariates.
			3. multivariate analysis of exposure on outcome, with the inclusion of covariates an confounders.
		2. one-sample MR
			1. estimation of d.hat (genotype predicted exposure)
			2. two-stage least square - multivariate analysis of d.hat (exposure) on outcome, with the inclusion of covariates.
			3. two-stage least square - multivariate analysis of d.hat (exposure) on outcome, with the inclusion of covariates and confounders. 
	
	arguments:
		1. workingdata: a data frame that includes all needed variables (grs, exposure, outcome, covariates, and confounders) in columns and individuals in rows.
		2. grs: a string that is the name of the column that holds the GRS (PRS) values.
		3. exposure: a string that is the name of the column that holds the exposure values.
		4. outcomes: a vector listing all of the column names one wishes to analyze as outcomes
		5. covariates: a vector listing all of the column names one wishes to analyze as covariates
		6. confounders: a vector listing all of the column names one wishes to analyze as confounders.
	
	example:
			my_1sam_mr = onesample_mr(
			 workingdata = mydata,
			 grs = "Yengo_Weighted_GRS", 
			 exposure = "bmi",
			 outcomes = traits,
			 covariates = c("visitdate","sex","age"),
			 confounders = c("PC3" , "smoking" , "diet" , "carbs_g") )

	returned: The function returns a data frame of sample size, effect estiamtes, standard error, and p-values. The following are column names in the data frame.
		Column name definitions:
			1. outcome: the outcome variable
			2. obs_uni_n: the sample size for the observation univariate analysis
			3. obs_uni_fstat: the observational univariate model Fstatistic
			4. obs_uni_adjRsq: the observational univariate total model adjusted R-squared
			5. obs_uni_VarExp: the observational univariate variance explained by exposure
			6. obs_univ_Beta: the observational univariate effect estimate
			7. obs_univ_SE: the observational univariate standard error
			8. obs_univ_Pval: : the observational univariate p-value
			9. obs_multi_n: the sample size for the observational multivariate analysis 
			10. obs_multi_fstat: the observational multivariate (w/covariates) Fstatistic
			11. obs_multi_adjRsq: the observational multivariate (w/covariates) total model adjusted R-squared
			12. obs_multi_VarExp: the observational multivariate (w/covariates) variance explained by the exposure
			13. obs_multi_Beta: the observational multivariate (w/covariates) effect estimate for the exposure
			14. obs_multi_SE: the observational multivariate (w/covariates) standard error for the exposure
			15. obs_multi_Pval: the observational multivariate (w/covariates) p-value for the exposure
			16. obs_multi_cfs_n: the sample size for the observational multivariate (w/covariates & confounders)
			17. obs_multi_cfs_fstat: the observational multivariate (w/covariates & confounders) Fstatistic
			18. obs_multi_cfs_adjRsq: the observational multivariate (w/covariates & confounders) total model adjusted R-squared
			19. obs_multi_cfs_VarExp: the observational multivariate (w/covariates & confounders) variance explained by the exposure
			20. obs_multi_cfs_Beta: the observational multivariate (w/covariates & confounders) effect estimate for the exposure
			21. obs_multi_cfs_SE: the observational multivariate (w/covariates & confounders) standard error for the exposure
			22. obs_multi_cfs_Pval: the observational multivariate (w/covariates & confounders) p-value for the exposure
			23. tsls2_n: the sample size for the two-stage least square (w/covariates) analysis
			24. tsls2_fstat: the two-stage least square (w/covariates) Fstatistic
			25. tsls2_adjRsq: the two-stage least square (w/covariates) total model adjusted R-squared
			26. tsls2_VarExp: the two-stage least square (w/covariates) variance explained by the exposure
			27. tsls2_Beta: the two-stage least square (w/covariates) effect estimate for the exposure on the outcome
			28. tsls2_SE: the two-stage least square (w/covariates) standard error for the exposure
			29. tsls2_Pval: the two-stage least square (w/covariates) p-value for the exposure
			30. tsls2_cf_n: the sample size for the two-stage least square (w/covariates & covariates) analysis
			31. tsls2_cf_fstat: the two-stage least square (w/covariates & covariates) Fstatistic
			32. tsls2_cf_adjRsq: the two-stage least square (w/covariates & covariates) total model adjusted R-square
			33. tsls2_cf_VarExp: the two-stage least square (w/covariates & covariates) variance explained by the exposure
			34. tsls2_cf_Beta: the two-stage least square (w/covariates & covariates) effect estimate for the exposure on the outcome
			35. tsls2_cf_SE: the two-stage least square (w/covariates & covariates) standard error for the exposure
			36. tsls2_cf_Pval: the two-stage least square (w/covariates & covariates) p-value for the exposure

### top_onesample_mr()
	name: top_onesample_mr()

	description: A funtion that takes the out put of the onesample_mr() function as input and returns all observations with a p-value smaller than or equal to alpha, in a simplified data frame. 

	arguments:
		1. mydataframe: the outout data frame from onesample_mr()
		2. alpha: the p-value threshold for screening observations

	example:
		top_outcomes = top_onesample_mr(my_1sam_mr, alpha = 0.05)

	returned: A simplified data frame with the two-stage least square results screened for. The data frame includes sample size, beta, se, p-value, and variance explained.
		1. N, sample size
		2. beta, TSLS effect estiamte
		3. se, TSLS standard error
		4. P, TSLS p-value
		5. eta-sq, variance explained for d.hat on outcome
	

### convert_2_longformat_onesample_mr()
	name: convert_2_longformat_onesample_mr()

	description: Converts the out data frame from the function onesample_mr() into a long format for ggplot.

	arguments: 
		1. data: the out data frame from the function onesample_mr()
		2. alpha: a p-value threshold to flag analysis as significant or not. Needed for the plotting function.

	example:
		MRplotdata = convert_2_longformat_onesample_mr(data = my_1sam_mr, alpha = BF)

	returned: A long formated data frame
		- Columns in data frame are:
			1. 	outcome
			2. n: sample size
			3. beta: effect estimate
			4. se: standard error
			5. pval: p-value
			6. analysis: analysis type from function onesample_mr()
	
### onesample_mr_forrest_plot()
	name: onesample_mr_forrest_plot()

	description: A function for making a forrest plot. The data input is the output from convert_2_longformat_onesample_mr().

	arguments:
		- data: the input is the output from convert_2_longformat_onesample_mr()
		- alpha: the p-value threshold for selection what to plot
		- analysis: the analysis type to screen p-values for
		- arrange_by: column names to sort the forrest plot by
		- title: A string that will be set as the title of the plot
		- covariates: a vector listing all of the covariates used in the function onesample_mr()
		- confounders: a vector listing all of the confounders used in the function onesample_mr()
		- outcome_label: a string that will become the outcome (y-axis) label
		- exposure_label: a string that will become the exposure (x-axis) label
		- filter_data = TRUE: FALSE if the data does not need to be filtered 'alpha'
		- column_4_x = "analysis": The column of data to set as the x-axis (this could be sex for example)
		- column_4_color = "analysis": The column of data to set as the color partitioning scheme

	example:
		P1 = onesample_mr_forrest_plot(data = MRplotdata,
				alpha = BF, analysis = "tsls",
				arrange_by = "pval",
				title = "Metabolite thats surpass BF correction (n = 39)",
				covariates = cvs,
				confounders = cnfs,
				outcome_label = "metabolites",
				exposure_label = "normalized standard deviation change per unit increase BMI" )

	returned: A ggplot forrest plot is returned. 
	