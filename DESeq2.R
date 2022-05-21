#TXimport DESeq2

#####################################################
#BiocManager::install("DESeq2")
#BiocManager::install("tximport")
#BiocManager::install("biomaRt")
#BiocManager::install("rhdf5")
#BiocManager::install("apeglm")
#####################################################

#TODO: SET TO YOUR WROKING DIRECTORY
setwd("/YOUR/FILE/PATH/HERE")
getwd()

###
### 0. Load Libraries
###
library(DESeq2)
library(tximport)
library(biomaRt)
library(rhdf5)
library(apeglm)
library(EnhancedVolcano)


###
### 1. Annotation
###

# Retrieve gene names from ENSEMBL using biomaRt:
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "mmusculus_gene_ensembl",
                         host = 'ensembl.org')

#TRANSCRIPT TO GENE
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "transcript_version", "ensembl_gene_id",
                                     "external_gene_name", "description", "transcript_biotype"), mart = mart)
head(t2g)

t2g <- dplyr::select(t2g, TXNAME = ensembl_transcript_id,
                     GENEID = ensembl_gene_id)

head(t2g)

###
### 2. Read Files
###

#Specify where the kallisto results are stored. 
sample_id <- dir(file.path("kallisto"))
sample_id

#A list of paths to the kallisto results indexed by the sample IDs is collated with:
kal_dirs <- file.path( 'kallisto', sample_id, 'abundance.h5')
kal_dirs

#The next step is to load an auxiliary table that describes the experimental design and the relationship between 
# the kallisto directories and the samples (called sample_to_condition):
s2c <- read.table(file.path( "/YOUR/FILE/PATH/HERE/metadata.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c

s2c <- dplyr::mutate(s2c, path = kal_dirs)
s2c

s2c$path

###
### 3. Tximport
###

txi_kallisto <- tximport(s2c$path, type = 'kallisto',
                         tx2gene = t2g, ignoreTxVersion = TRUE) #, txOut = TRUE
head(txi_kallisto$counts)

#Get TPM
tpm <- txi_kallisto$abundance
colnames(tpm) <- s2c$sample
dim(tpm)

tpm <- as.data.frame(tpm)

#TODO: SET TO YOUR CONDITIONS
control<- which(s2c$condition=="BV2_HA" )
control
overexpressed<- which(s2c$condition=="BV2_Rhox5")
overexpressed

#Take mean across replicates
tpm$ctrl_mean <- rowMeans(tpm[ , c(control)], na.rm=TRUE)
head(tpm )
tpm$rhox_mean <- rowMeans(tpm[ , c(overexpressed)], na.rm=TRUE)
head(tpm )
#convert rownames to first column
library(dplyr)
tpm <- tibble::rownames_to_column(tpm, "target_id")

#Prepare for DESeq2
dds_kallisto <- DESeqDataSetFromTximport(txi_kallisto, s2c, ~condition)

#Filter for lowly expressed genes
threshold <- 10
keep <- rowSums(counts(dds_kallisto)) >= threshold  
dds_kallisto <- dds_kallisto[keep, ]
print(paste('dimensions of filtered genes', dim(dds_kallisto)[1], dim(dds_kallisto)[2]))

###
### 4. Differential Expression Analysis
###

dds_kallisto <- DESeq(dds_kallisto)

resultsNames(dds_kallisto) # lists the coefficients

# Results tables are generated using the function results, which extracts a results table with log2 fold changes, 
# p values and adjusted p values. 

#By specifying the lfcThreshold, we incorporate the fold change criterion into the Wald statistic, and are now testing our observed result against the values we would 
# expect under the null distribution where the difference between the two groups is no larger than the lfc specified. This is a much more conservative threshold, 
# and you should expect far fewer genes. 
res <- results(dds_kallisto, name="SET TO YOUR COEFFICIENT")
res

summary(res)

#Further Filtering
#filt1 = res[which(res$pvalue < 0.05), ]
#filt2 = filt1[which(filt1$padj < 0.1), ]
#print(paste('DEGs found', dim(filt2)[1], sep=' '))

#saveRDS(res, file = 'DESeq2_BV2_Rhox5.rds')

#Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes
resLFC <- lfcShrink(dds_kallisto, coef="SET TO YOUR COEFFICIENT", type="apeglm")
resLFC

#Order the results table by smallest p-value and view summary:
resOrdered <- res[order(res$pvalue),]
summary(res)
#how many with p<0.01?
sum(res$padj < 0.1, na.rm=TRUE)

# transformation of count data for visualization
# blind is whether the transformation should be blind to the sample information specified by the design formula
vsd <- vst(dds_kallisto, blind=FALSE)
head(assay(vsd), 3)

#Heatmap
dev.off()
library(ggplot2)
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition,'-', vsd$sample)
colnames(sampleDistMatrix) <- paste(vsd$condition,'-', vsd$sample)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
library(pheatmap)
png(filename = "heatmap.png",
    width = 190, height = 150, units = "mm", res=400)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

#PCA
plotPCA(vsd, "condition")
ggsave("pca.png", width = 190, height = 150, units ="mm")

#Volcano
dev.off()
png(filename = "volcano.png",
    width = 190, height = 150, units = "mm", res=400)
EnhancedVolcano(resLFC,
                title = NULL,
                caption = NULL,
                subtitle = NULL,
                legendPosition = 'bottom',
                lab = rownames(resLFC),
                x = 'log2FoldChange',
                y = 'pvalue',
                labSize= 0,
                pCutoff = 0.01, 
                FCcutoff =1)
dev.off()


#Likelihood Ratio Test
# A simple likelihood ratio test where the full design was ~condition:
dds_lrt <- DESeq(dds_kallisto, test="LRT", reduced=~1)
res_lrt <- results(dds_lrt)
res_lrt

###
### 5. Write Results to File
###

#Wald
#if I filter, use res here
DESeq2_kallisto <- results(dds_kallisto)
#summary(DESeq2_kallisto)
DESeq2_kallisto <- as.data.frame(DESeq2_kallisto)
DESeq2_kallisto <- dplyr::mutate(DESeq2_kallisto, target_id = rownames(DESeq2_kallisto))

#alpha <- 0.05
#significant_kallisto <- dplyr::filter(DESeq2_kallisto, padj < alpha)
#write.table(DESeq2_kallisto,file= "BV2_Rhox5_DESeq2.tsv", quote=FALSE, sep='\t', row.names = FALSE )

#LRT
res_lrt <- as.data.frame(res_lrt)
res_lrt <- dplyr::mutate(res_lrt, target_id = rownames(res_lrt))
#significant_lrt <- dplyr::filter(res_lrt, padj < alpha)
#write.table(res_lrt,file= "BV2_Rhox5_DESeq2_LRT.tsv", quote=FALSE, sep='\t', row.names = FALSE )

#Merge Wald and LRT Tables
results_table <- merge(DESeq2_kallisto, res_lrt, by.x  = 'target_id',by.y = 'target_id', sluth_tablert = FALSE)
head(results_table)

#Merge this with TPM
merged <-merge(results_table,tpm,by=c("target_id"))
write.table(merged,file= "DESeq2_MERGE.tsv", quote=FALSE, sep='\t', row.names = FALSE )
sessionInfo()

