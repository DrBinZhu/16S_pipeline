# Automatic analysis of 16S rRNA sequencing data edited by Dr. Bin Zhu at VCU.

# How to use
Please see the example files and modify the formate of your input, change the names of the input files to 'reads_table.csv' and 'metadata.csv'. Some Windows system require moving the '16s_analysis_pipeline_BZ' folder out of download folder before running the scripts. Open '16s_analysis_pipeline_*.R' by Rstudio, select all the commands, and click run.

Sample_IDs in 'reads_table.csv' can be put in either rows or columns, and the order of the Sample_IDs are not necessary to be the same as those in 'metadata.csv'. Additional samples in 'reads_table.csv' do not impact the analysis. Just make sure that all the samples in 'metadata.csv' exist in 'reads_table.csv'.

Several parameters can be changed in the 'parameter_list.csv' to modify your output.

# Methods
Taxa with at least 0.1% (or 0.01%) relative abundance in at least 5% (or 15%), respectively, of the samples are kept in the feature table of the 16S rRNA profiles. Samples in the feature table with a total sample reads less than 5,000 are excluded. 

The feature table of the 16S rRNA profiles is normalized by rarefaction to the depth of the lowest number of reads in a sample (>5,000). Alpha diversity is quantified by calculating Shannon index, evenness, and observed taxa number using the ‘vegan’ package in R. The difference in the alpha diversity is measured by the Mann-Whitely U test. Beta diversity is measured and visualized by a Non-metric Multidimensional Scaling (NMDS) of Bray-Curtis distances using the ‘vegan’ package in R. The difference in the beta diversity is tested by the ‘adonis2’ function in the ‘vegan’ package in R. 

The differential abundance ananlysis is performed using the 'ALDEx2' package with P-values tested by the Mann-Whitely U test and adjusted by the Benjamini-Hochberg Procedure. The fold changes and adj-P-values are shown in the columns 'diff.btw' and 'wi.eBH', respectively, in the output file 'abundance_difference.csv'.

To generate a heatmap, the feature table of the 16S rRNA profiles is converted to a relative abundance table. Samples in the relative abundance table are clustered by calculating the Euclidean distance and using the ‘complete’ method in the ‘pheatmap’ function in R. 

To measure bacterial association, the feature table of the 16S rRNA profiles is normalized by the centered log-ratio transformation. The association between each taxa dyad is measured by the Spearman’s correlation. The coefficient and significance of a correlation are quantified by an R- and P-value, respectively. The P-values are adjusted by the Benjamini and Hochberg procedure. The R-values with matched adjusted P-values larger than 0.05 are replaced by zeros to remove insignificant correlations in the next step. Taxa in the VMB are clustered according to the adjusted R-values of the Spearman’s correlation using the ‘pheatmap’ function with default settings in R. For a network with only abundant taxa, only taxa with at least 10% (or 1%) relative abundance in at least 5% (or 15%), respectively, of the samples are kept in the feature table of the 16S rRNA profiles. All other steps are the same.

# Cite this paper
Zhu, B., Diachok, C., Edupuganti, L., Edwards, D. J., Donowitz, J. R., Tossas, K., ... & Buck, G. A. The Utility of Voided Urine Samples as a Proxy for the Vaginal Microbiome and for the Prediction of Bacterial Vaginosis. Infectious Microbes & Diseases, 10-1097.
