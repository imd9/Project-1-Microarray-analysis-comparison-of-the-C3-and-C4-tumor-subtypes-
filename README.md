# Project Description
Colorectal cancer (CRC) is the third most common type of cancer and the fourth most common cause of death in the world. Pathological staining is the most widely used method to determine the presence of CC, however, it does not accurately predict recurrence of CRC. About 20% of patients who are diagnosed with stage II or III CRC develop recurrence. Researchers have looked into Gene Expression Profiles (GEPs) through the use of microarrays. The Marisa et al study established a classification of the CC subtypes based on their molecular features by exploiting “genome-wide mRNA expression analysis” through the use of microarrays. Initially only three subtypes of CC were identified, but through the Marisa et al study, six subtypes were classified, more accurately reflecting the molecular heterogeneity of CC. 

This analysis focuses only on reproducing the results from the comparison of the C3 and C4 tumor subtypes from Marisa et al. 134 samples were analyzed and taken as an input. 

Referenece: Marisa et al. Gene Expression Classification of Colon Cancer into Molecular Subtypes: Characterization, Validation, and Prognostic Value. PLoS Medicine, May 2013. PMID: 23700391

# Contributors

Monica Roberts: Data curator (monicapr@bu.edu)

Preshita Dave: Programmer (preshita@bu.edu)

Italo Duran: Analyst (duran01@bu.edu)

# Repository Contents
# Usage 
R script.R

# Programmer
The programmer.R script executes the necessary commands for assessing the quality of the data as well as preprocessing steps such as normalization of the data through the RMA method. 

Input: 

134 CEL files

Output: 

RMA Normalized and batch effect corrected expression intensity values

# Analyst
The analysis_final.R script performs noise filtering and dimension reduction (based off on the PCA results in the programmer.R script), followed by hierarchial clustering and subtype discovery. 

Input: 

Output after running the programmer.R script

Output: 

Gene matrices after passing the data through 3 filters

Expression matrices after passing various filters from the chi-square and the Welch-t test

Heatmap to visualize the differentially expressed genes in the C3 and C4 subtypes


