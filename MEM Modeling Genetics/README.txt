#R SESSION INFORMATION#
R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
[4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggplot2_3.4.0         reshape_0.8.9         RColorBrewer_1.1-3    matrixStats_0.63.0    plyr_1.8.8           
 [6] viridis_0.6.2         viridisLite_0.4.1     circlize_0.4.15       ComplexHeatmap_2.12.1 emmeans_1.8.3        
[11] lmerTest_3.1-3        lme4_1.1-31           Matrix_1.5-3   

#RUNNING GENOMIC AND IMAGING MIXED EFFECT MODELING (MEM)#
Files labeled "...ANOVAfinal.R" are used for mixed effect modeling with categorical variables including genomics and transcriptomic pathway class. Run each script within its own subfolder to prevent outputs being overwritten.
Run times are on the order of minutes

Inputs for these scripts are:
MultiReg_imaging_genetics.csv
-includes genomic data and imaging values for samples with this data 
combinations.csv
-lists every combination of imaging variable and contrast enhancement with the transformation used to meet the assumption of equal variance
The transcriptomic pathway script also requires:
"Garofanoclass.csv"
-lists the transcriptomic class of all samples classified 

Each script follows the following protocol: 
1) read in dataand libraries, preprocess/ labeling, and creation of empty variables to store results 
2)The for loop iterates over the combinations specified in the combinations file which says which imaging variable is being modeled and if NE, CE, or all (Total) samples are being examined
3)for the first run of the script all combinations are modeled with no variable transformation ("None"). After running, the quality control plots "....QC.png" are examined to ensure the assumption of equal variance is 
met which is required for ANOVA. If the plots show variance bias, the "combinations.csv" file is edited for the relevant combination to "Log10", "Sqrt", or "Recip" until the variance reasonably equalized for each combination.
*final transforms used for all of our models are in the "combination_transforms" folder. 
4)MEM model is constructed and anova analysis is run on said model. Pertinent outputs from this analysis are stored in previously created empty variables.
5)A heatmap of all p-values for pairwise comparison between geneotypes is saved as "...pval_ht.png"
6) A smaller heatmap "...pval_ht_subset.png" that trims off rows with values very close to 1 (not significant) for better visibility of significant pariwise comparisons
7) genotypes that had significant pairwise results (pvalue<0.05) were analyzed to see if they showed significance against all other genotypes in a nested for loop
8) values of this significant genotype versus the rest of the samples were saved as "....specifics.png"
9) aggregate files are saved for comparing across imaging values and combinations 

Outputs for each combination (loop iteration):
"..._anova.csv": a summary of results from the MEM model and subsequent ANOVA
"...._constrasts.csv": p-values for contrasts between each examined gene in the MEM
"...pval_ht_subset.png" and "...pval_ht.png": heatmaps of p-values of pairwise comparisons between genotypes
"...pvalues.png": p-values of pairwise comparisons between genotypes
"....QC.png": plot for examining if the model meets the assumption of equal variance 

Outputs if a significant genotypes is found for a combination:
"_specifics.png": comparison of MEM corrected means for significant genotypes against all other genotypes

Summary Outputs:
...._percentvariance_attributed.csv: percent variance attributed to each fixed genetic term within the MEM for each imaging variable in all three contexts (CE,NE,Total)
...._pvalue.csv: level of significance for each fixed genetic term within the MEM for each imaging variable in all three contexts (CE,NE,Total)