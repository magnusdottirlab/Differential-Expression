# Import results from Sleuth and DESeq2 and make Rank files
library(tidyverse)

###
### 0. Import Datasets
###

sleuth <- read.delim(file.path("/YOUR/FILE/PATH/HERE/sleuth_table_gene.tsv"), header = TRUE, sep = "\t")    #, fill = TRUE)
deseq <- read.delim(file.path("/OUR/FILE/PATH/HERE/DESEQ.tsv"), header = TRUE, sep = "\t")

#keep only needed columns 
sleuth <- sleuth[c(1,2,6,7)]
deseq <- deseq[c(1,3,6,7,9,12,13)]

###
### 1. Compare Results
###

# Compare the results of both tools

# Merge these into one table, keeping only the common genes:
DEG <- merge(sleuth,deseq,by.x="target_id",by.y="target_id",all =FALSE)
setwd("/YOUR/FILE/PATH/HERE")
getwd()
#write.table(DEG,file= "DEG.tsv", quote=FALSE, sep='\t', row.names = FALSE )


#Now we can determine the number of shared and unique significant genes:
print("Number of Differentially Expressed Genes Detected by DESeq2:")
nrow(deseq)
print("Number of Differentially Expressed Genes Detected by Sleuth:")
nrow(sleuth) 
print("Number of Differentially Expressed Genes Detected by Both:")
nrow(DEG)

#number of common elements
length(intersect(sleuth$target_id, deseq$target_id))

# Setdiff indicates which elements of a vector or data frame X are not existent in a vector or data frame Y
#in sleuth not deseq
in_sleuth_not_deseq <- setdiff(sleuth$target_id, deseq$target_id)
length(in_sleuth_not_deseq)
in_sleuth_not_deseq <- as.data.frame(in_sleuth_not_deseq)

#in deseq not sleuth
in_deseq_not_sleuth <- setdiff(deseq$target_id, sleuth$target_id)
length(in_deseq_not_sleuth)
in_deseq_not_sleuth <- as.data.frame(in_deseq_not_sleuth)

###
### 2. Calculate rank statistic AGGREGATED
###

# Rank = sign(log2(fc))*-log10(qValue)

#https://www.marsja.se/r-add-column-to-dataframe-based-on-other-columns-conditions-dplyr/
DEG_AGG <- DEG %>% rowwise() %>%
  mutate(RANK_AGG = sign(log2FoldChange.x)*(-log10(qval)))
nrow(DEG_AGG)

sapply(DEG_AGG, class) 

#List by rank (Descending)
DEG_AGG <- DEG_AGG[
  with(DEG_AGG, order(-RANK_AGG)),
]

#Delete rows with empty ext_gene name
DEG_AGG <- as.data.frame(DEG_AGG)
DEG2_AGG<-DEG_AGG[!(DEG_AGG$ext_gene==""),]

###
### 3. Write to file AGGREAGTED
###

#https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#RNK:_Ranked_list_file_format_.28.2A.rnk.29

#Only want gene and rank
DEG3_AGG <- DEG2_AGG[c(2,11)]

#Save to respective folder in GSEA 
getwd()
setwd("/YOUR/FILE/PATH/HERE/")
getwd()

write.table(DEG3_AGG,file= "AGGREGATED.rnk.txt", quote=FALSE, sep='\t', row.names = FALSE, col.names = FALSE )

