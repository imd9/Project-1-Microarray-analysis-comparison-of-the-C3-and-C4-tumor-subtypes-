# Part 4. Noise Filtering & Dimensionality Reduction

library(gplots)
# making a data-frame from loading data in:
data<-data.frame(read.csv('/projectnb/bf528/users/swiss_cheese_2022/project_1/project-1-swisscheese/corrected_batch_values.csv', sep=' ', header=TRUE))

#1. Expressed in at least 20% of samples (i.e. for each gene, at least 20% of the gene-expression values must be > log2(15)
expressed_20<-data[rowSums(data > log2(15)) >= (0.2*ncol(data)), ]

# 2. Have a variance significantly different from the median variance of all probe sets using a threshold of p < 0.01:
# Two-tailed 
deg_free<-ncol(expressed_20)-1
chi_low<-qchisq(0.01/2, deg_free)
chi_up<-qchisq(1-0.01/2, deg_free)

# The test statistic for each gene:
expressed_20$variance <- apply(expressed_20, 1, var)
expressed_20$test_stat <- deg_free*expressed_20$variance/median(expressed_20$variance)
chi_sq<-subset(expressed_20, test_stat > chi_up | test_stat < chi_low)
# Removing extra columns:
chi_sq<-subset(chi_sq, select = -c(variance, test_stat))

# 3. Have a coefficient of variation > 0.186 :
variation<-subset(chi_sq, apply(chi_sq, 1, function(x) sd(x)/mean(x)) > 0.186)

# Deliverables:

# 4. Write out a different file containing the gene expression matrix for genes passing all three of the filters from 4.1, 4.2, and 4.3 :
write.csv(variation, '/projectnb/bf528/users/swiss_cheese_2022/project_1/project-1-swisscheese/all_filtered_expression_matrix.csv')
# 5) Write out the expression matrix for probe sets that pass the expression threshold(biologist part)
write.csv(chi_sq, '/projectnb/bf528/users/swiss_cheese_2022/project_1/project-1-swisscheese/differential_filtered_expression_from_matrix.csv')
# Report the number of genes that pass all of these thresholds:
print(paste0('The number of genes that pass all of these thresholds: ', nrow(variation)))


##########################################################

# Part 5. Hierarchical Clustering & Subtype Discovery


# 1. Perform hierarchical clustering on your fully filtered data matrix from Part 4.4.
  # Be sure to check that you are clustering the patients and not the genes.
welch_t_test<-function(data_frame, file_name) {
  dendrogram<-hclust(dist(t(data_frame)))
  
# 2. Cut the dendrogram such that the samples are divided into 2 clusters. How many samples are in each cluster?
  clusters<-cutree(dendrogram, 2)
  print(paste0('Number of samples in cluster 1: ', sum(clusters==1)))
  print(paste0('Number of samples in cluster 2: ', sum(clusters==2)))
  
# 4. Using the expression matrix from Part 4.4 and the cluster memberships from Part 5.2,
# identify genes differentially expressed between the two clusters using a Welch t-test (results in a ranked list).
# Write out a data-frame containing the probeset ID, t-statistic, p-value, 
# and adjusted p-value (i.e. FDR, see the p.adjust function) to a comma separated file for each comparison.
# How many genes are differentially expressed at adjusted 𝑝<0.05 between the clusters for both lists?
  
  cluster1<-data_frame[, clusters==1]
  cluster2<-data_frame[, clusters==2]
  
  # Empty list to fill with data-frames for all the genes.
  t_test<-vector('list', nrow(data_frame))
  for (row in 1:nrow(data_frame)) {
    # Welch t-test:
    welch_t_test2<-t.test(cluster1[row,], cluster2[row,])
    # Write out a data-frame containing the probeset ID, t-statistic, p-value, and adjusted p-value.
    t_test[[row]]<-data.frame(t=welch_t_test2$statistic, p=welch_t_test2$p.value)}
  # Putting it into one data-frame:
  t_test<-do.call(rbind, t_test)
  # Setting gene names as row names:
  row.names(t_test)<-row.names(data_frame)
  # adjusted p value:
  t_test$p_adj<-p.adjust(t_test$p, method='fdr')
  # Write out a data-frame comma separated file:
  write.csv(t_test, file_name)
  
  # Print statements
  pval_adj<-subset(t_test, p_adj < 0.05)
  print(paste0('Number of differentially expressed genes at 𝑝<0.05 between the two clusters: ', nrow(pval_adj)))
  print('The genes that were representative of each cluster were determined by t-test.')
  print('The top genes with the largest differences in the clusters are: ')
  print(paste0(head(rownames(t_test[rev(order(abs(t_test$t))),]))))
}

# Perform the t-test analysis (Biology part)
print('T-test for expression matrix with all filters')
welch_t_test(variation, '/projectnb/bf528/users/swiss_cheese_2022/project_1/project-1-swisscheese/all_filters_t_test.csv')
print('T-test for expression matrix for probesets that pass the expression threshold')
welch_t_test(chi_sq, '/projectnb/bf528/users/swiss_cheese_2022/project_1/project-1-swisscheese/differential_filter_t_test.csv')

# 3. Create a heatmap of the gene-expression of each gene across all samples.
# Load metadata
metadata <- read.csv('/project/bf528/project_1/doc/proj_metadata.csv')
metadata <- subset(metadata, select = c(geo_accession, cit.coloncancermolecularsubtype))
#metadata <- metadata[metadata$geo_accession %in% colnames(variation),]
colors <- ifelse(metadata$cit.coloncancermolecularsubtype == 'C3', 'red', 'blue')

# Save heatmap as a file
png('/projectnb/bf528/users/swiss_cheese_2022/project_1/project-1-swisscheese/heatmap.png', width=1920, height=1080, res=100)
heatmap.2(as.matrix(variation), xlab='Patient Samples', ylab='Genes', 
          main='Gene Expression',
          ColSideColors = colors, trace='none', density.info = 'none',
          key.xlab='Gene Expression Level', scale='row', margins=c(15,7))

# Legend for Molecular Subtype
legend(x=0.9, y=1, xpd=TRUE, inset=c(-0.15,0),
       legend=c('C3', 'Other'), title='Molecular Subtype', fill=c('red', 'blue'))
dev.off()
