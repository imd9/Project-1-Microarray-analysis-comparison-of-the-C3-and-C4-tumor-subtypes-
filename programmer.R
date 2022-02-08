# loading installed packages
library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)

# loading the data by setting the working directory
setwd('/usr4/bf527/preshita/project-1-swisscheese/samples')
celpath = getwd()
setwd('/usr4/bf527/preshita/project-1-swisscheese')

# reading the data 
data = ReadAffy(celfile.path = celpath)

# 
data.rma = rma(data)

# creating an expression dataset which can be written out as a csv 
data.matrix = exprs(data.rma)
# write.table(data.matrix, file="rma.txt", quote=FALSE,sep="\t")

# Pset takes too much time

# https://www.bioconductor.org/packages/devel/bioc/vignettes/affyPLM/inst/doc/QualityAssess.pdf
Pset <- fitPLM(data, normalize=TRUE, background=TRUE)
# RLE(Pset,main="RLE for CRC dataset")


stats_RLE<-data.frame(RLE(Pset,type="stats"))
stats_RLE<-as.data.frame(t(stats_RLE))
jpeg("hist_RLE.jpg")
hist(stats_RLE$median, xlab="median RLE values", main = "Histogram plot of median RLE values", col='lightblue')
dev.off()


stats_NUSE<-data.frame(NUSE(Pset,type="stats"))
stats_NUSE<-as.data.frame(t(stats_NUSE))
jpeg("hist_NUSE.jpg")
hist(stats_NUSE$median, xlab="median NUSE values", main = "Histogram plot of median NUSE values", col = 'lightgreen')
dev.off()

# Combat for batch effect correction 
metadata=read.csv('/project/bf528/project_1/doc/proj_metadata.csv')
batch_effect_correction <- ComBat(dat = data.matrix, batch = metadata$normalizationcombatbatch, mod = model.matrix(~as.factor(metadata$normalizationcombatmod), data=metadata))
write.table(batch_effect_correction, file="/usr4/bf527/preshita/project-1-swisscheese/corrected_batch_values.csv")


#PCA 




