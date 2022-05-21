#Sleuth Gene-Level Analysis

library("sleuth")
packageVersion("sleuth")
library("dplyr")
library('cowplot')
library('shiny')
library('ggplot2')

#TODO: SET TO YOUR OWN WORKING DIRECTORY
setwd("FILE/PATH/HERE")
getwd()

#Specify where the kallisto results are stored. 
sample_id <- dir(file.path("kallisto"))
sample_id

kal_dirs <- file.path( "kallisto", sample_id)
kal_dirs

#TODO: SET TO YOUR OWN FILE PATH
s2c <- read.table(file.path( "FILE/PATH/HERE/metadata.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c

s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)

#Collect gene names with:
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "mmusculus_gene_ensembl",
                         host = 'ensembl.org')

t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name", "description"), mart = mart)


t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
head(t2g)

so_rhox <- sleuth_prep(s2c, target_mapping = t2g, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

so_rhox <- sleuth_fit(so_rhox, ~condition, 'full')
so_rhox <- sleuth_fit(so_rhox, ~1, 'reduced')
so_rhox <- sleuth_lrt(so_rhox, 'reduced', 'full')

#Obtaining gene-level differential expression results
# When running the command ‘sleuth_results,’ sleuth uses the p-values from comparing transcripts to make a gene-level determination and
# perform gene differential expression.

sleuth_table_gene <- sleuth_results(so_rhox, 'reduced:full', 'lrt', show_all = FALSE)
#sleuth_table_gene <- dplyr::filter(sleuth_table_gene, qval <= 0.05)
#head(sleuth_table_gene, 20)

write.table(sleuth_table_gene,file= "sleuth_table_gene.tsv", quote=FALSE, sep='\t', row.names = FALSE )

saveRDS(so_rhox, file = 'so_GENE.rds')

sleuth_live(so_rhox)

#####################
# Manually download Summaries -> Processed Data

### Plotting
## Maps
# PCA
plot_pca(so_rhox, color_by = 'condition', units = "est_counts")  
ggsave("pca_est_gene.png", width = 190, height = 150, units ="mm")
plot_pca(so_rhox, color_by = 'condition', units = "tpm")  
ggsave("pca_tpm_gene.png", width = 190, height = 150, units ="mm")

# Heatmap
png(filename = "heatmap_gene.png",
    width = 190, height = 150, units = "mm", res=400)
plot_sample_heatmap(so_rhox)
dev.off()

