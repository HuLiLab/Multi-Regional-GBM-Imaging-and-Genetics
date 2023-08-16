Author: Taylor Weiskittel
Contact: tweiskittel94@gmail.com
R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
[4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] reshape2_1.4.4 MuMIn_1.47.1   lmerTest_3.1-3 lme4_1.1-31    Matrix_1.5-3  

#RUNNING PATHWAY SIGNATURES AND IMAGING MIXED EFFECT MODELING (MEM)#
"imaging_pathclass_tables.R" is a script that creates tables of summary statistics and MEM model parameters for each imaging parameter and biological signature pair in the CE,NE, and both (Total) 

Inputs:
MultiReg_advimaging_pathclass.csv and MultiReg_imaging_pathclass.csv: contains imaging data and pathway signature scores for each sample

Protocol: 
1) read in data and libraries
2) create new matrices with CE only and NE only samples 
3) scale all imaging parameters and biological signatures for the Total, NE, and CE data matrices 
4) loop through each pathway classification type 
5) nested loop through each imaging variable 
6) for each iamging variable calculate the correlation and correlation p-value for CE, NE, and Total sample sets. Then create a model for the three sample sets that models the biologicall signature as a fixed effect and patients as 
random effects for the prediction of each imaging variable. Record the R squared, coefficient associated with the biological signature, and SE for this coefficient. 
7) save the result in pathway classification specific tables 
8) repeat 1-7 for the second set of imaging features ("MultiReg_imaging_pathclass.csv")

Outputs for each Pathway class and imaging feature set:
"...._statistics.csv": results matrix of the statistical and MEM model measures for a garofano class  


