# Automatic analysis on 16s rRNA sequencing data edited by Dr. Bin Zhu at VCU.

# How to use
Please see the example files and modify the formate of your input, change the names to the input files to 'reads_table.csv' and 'metadata.csv'. Windows users need to move the '16s_analysis_pipeline_BZ' folder out of download folder before running the scripts. Open '16s_analysis_pipeline_*.R' by Rstudio, select all the commands, and click run. You need to input a control name for each comparison.

# Methods
Taxa with at least 0.1% (or 0.01%) relative abundance in at least 5% (or 15%), respectively, of the samples will be kept in the feature table of the 16S rRNA profiles. Samples in the feature table with a total sample reads less than 5,000 will be excluded. 

The feature table of the 16S rRNA profiles will be normalized by rarefaction to the depth of the lowest number of reads in a sample (>5,000). Alpha diversity will be quantified by calculating Shannon index, evenness, and observed taxa number using the ‘vegan’ package in R. The difference in the alpha diversity will be measured by the Wilcoxon test. Beta diversity will be measured and visualized by a Non-metric Multidimensional Scaling (NMDS) of Bray-Curtis distances using the ‘vegan’ package in R. The difference in the beta diversity will be tested by the ‘adonis2’ function in the ‘vegan’ package in R. 

To generate a heatmap, the feature table of the 16S rRNA profiles will be converted to a relative abundance table. Samples in the relative abundance table will be clustered by calculating the Euclidean distance and using the ‘complete’ method in the ‘pheatmap’ function in R. 

To measure bacterial interaction, the feature table of the 16S rRNA profiles will be normalized by the centered log-ratio transformation. The association between each taxa dyad will be measured by the Spearman’s correlation. The coefficient and significance of a correlation will be quantified by an R- and P-value, respectively. The P-values will be adjusted by the Benjamini and Hochberg procedure using the ‘adjust.p’ function in the ‘cp4p’ package in R. The R-values with matched adjusted P-values larger than 0.05 will be replaced by zeros. Taxa in the VMB will be clustered according to the adjusted R-values of the Spearman’s correlation using the ‘pheatmap’ function with default settings in R. For a network with only abundant taxa, only Taxa with at least 10% (or 1%) relative abundance in at least 5% (or 15%), respectively, of the samples will be kept in the feature table of the 16S rRNA profiles. All other steps are the same. 
