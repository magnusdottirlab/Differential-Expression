#Sleuth Transcript-Level Analysis

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

#A list of paths to the kallisto results indexed by the sample IDs is collated with:
kal_dirs <- file.path( "kallisto", sample_id)
kal_dirs

#The next step is to load an auxiliary table that describes the experimental design and the relationship between 
# the kallisto directories and the samples:
#TODO: SET TO YOUR OWN PATH
s2c <- read.table(file.path( "FILE/PATH/HERE/metadata.txt"), header = TRUE, stringsAsFactors=FALSE)

# By default, R converts any string variables (text, rather than numbers) into a special data type called factors. A factor is R’s way 
# of representing categorical variables so that they can be used in statistical models. However, our sample and path columns are not 
# categorical variables like this, so we tell R not to convert strings to factors with stringsAsFactors=FALSE. 

s2c

#Now the directories must be appended in a new column to the table describing the experiment. This column must be labeled 
#path, otherwise sleuth will report an error. This is to ensure that samples can be associated with kallisto quantifications.
s2c <- dplyr::mutate(s2c, path = kal_dirs)

#It is important to check that the pairings are correct:
print(s2c)

#Including gene names into transcript-level analysis
# In reading the kallisto output sleuth has no information about the genes that transcripts are associated with, but this can 
# be added allowing for searching and analysis of significantly differential transcripts by their associated gene names.

# We will add gene names from ENSEMBL using biomaRt (there are other ways to do this as well):

#Collect gene names with:
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "mmusculus_gene_ensembl",
                         host = 'ensembl.org')

t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name", "description"), mart = mart)
#"transcript_version" 

t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
head(t2g)

# Next, the "sleuth object" can be constructed. This object will store not only the information about the experiment, but also 
# details of the model to be used for differential testing, and the results. It is prepared and used with four commands that: 
# (1) load the kallisto processed data into the object 
# (2) estimate parameters for the sleuth response error measurement (full) model 
# (3) estimate parameters for the sleuth reduced model, and 
# (4) perform differential analysis (testing) using the likelihood ratio test. 

# This is where the model will be fit. We fit data from different conditions to distributions and then compare those distributions. Normal 
# distributions don't fit RNA-seq data well as raw read counts are discrete and cannot be below zero and a normal distribution is continuous with 
# negative values. Therefore, in order to use normal distribution techniques, some tools (including Sleuth) choose to 
# convert RNA-Seq data into a distribution that approximates the normal distribution. How can we approximate a normal distribution from negative 
# binomially distributed data? We can take logarithms of our abundance values, which makes them roughly normal. Note that while the log conversion 
# does reduce the long positive tail that the negative binomial distribution usually has, it can produce a skew to the left, especially at low mean
# and size values. It is not a perfect transformation to the normal by any means, but it has proven to be close enough in practice for most genes.

# The sleuth object must first be initialized with:
so_rhox <- sleuth_prep(s2c, target_mapping = t2g, extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

# We used sleuth_fit to fit two models, one including condition (called full) and one without (called reduced)
# The aim of fitting these two models is to test for an effect of condition; we want to assess whether the full model fits the data better 
# than the reduced model. If it is only a slightly better fit, condition probably doesn’t have much effect; if the model fits a lot better, 
# condition is probably important.

#  Then the full model is fit with:
so_rhox <- sleuth_fit(so_rhox, ~condition, 'full')

# What this has accomplished is to “smooth” the raw kallisto abundance estimates for each sample using a linear model with a parameter that
# represents the experimental condition. To test for transcripts that are differential expressed between the conditions, sleuth performs a 
# second fit to a “reduced” model that presumes abundances are equal in the two conditions. To identify differential expressed transcripts
# sleuth will then identify transcripts with a significantly better fit with the “full” model.

# The "reduced" model is fit with:
so_rhox <- sleuth_fit(so_rhox, ~1, 'reduced')

# We can test this with a likelihood ratio test, run using sleuth_lrt. It takes a Sleuth object (so), a null model (reduced in our case) 
# and an alternative model (full here). It then calculates how many times more likely the alternative model is than the null model 
# (the likelihood ratio), and calculates a p value, which is the probability of seeing this likelihood ratio if the variables added 
# in the full model have no effect.
so_rhox <- sleuth_lrt(so_rhox, 'reduced', 'full')
# To see the results of this test, open the Sleuth web interface with sleuth_live(so). Sleuth shows Wald test results by default, 
# so first go to the settings menu item and change test type to ‘likelihood ratio’. Now go back to analyses → test table. This table 
# now shows likelihood ratio test results, indicating which genes have a significantly better model fit when condition is included, 
# and so are likely to be affected by condition.

#The models that have been fit can always be examined with the models() function.
#TODO: USE THIS OUTPUT AS INPUT TO OTHER FUNCTIONS LATER. Mine is called 
models(so_rhox)

# All we have done so far is tested if our condition  has an effect on some of our genes. We have not yet quantified the effect size. 
# We can use a Wald test to measure the effect of a particular condition against the control condition for a variable. In the command below, 
# we test for an effect of condition, that is, the condition overexpressed or control. Wald tests also calculate p
# values to estimate the significance of any difference between conditions. Wald test p values represent the probability of seeing this data 
# if the true effect size for the condition was zero. So if the p value is very small, it is likely that the condition has a non-zero effect.

# We run a pairwise Wald test to compare the two different conditions with the sleuth_wt function. This test produces a
# list of genes that are deferentially expressed between the wild type and mutant samples:

so_rhox <- sleuth_wt(so_rhox, which_beta= "CHANGE_TO_YOUR_MODEL")
sluth_resultsWT <- sleuth_results(so_rhox, 'CHANGE_TO_YOUR_MODEL', show_all = FALSE)

# Calculate effect size based on log2 instead of natural log to make it more interpretable. 
sluth_resultsWT$raw_b <- exp(sluth_resultsWT$b)
sluth_resultsWT$log2_b <- log2(sluth_resultsWT$raw_b)
sluth_resultsWT$raw_se_b <- exp(sluth_resultsWT$se_b)
sluth_resultsWT$raw_log2_se_b <- log2(sluth_resultsWT$raw_se_b)

sleuth_table_tx <- sleuth_results(so_rhox, 'reduced:full', 'lrt', show_all = FALSE)

results_table <- merge(sleuth_table_tx, sluth_resultsWT, by.x  = 'target_id',by.y = 'target_id', sluth_tablert = FALSE)
head(results_table)

#--------------- Write TPM to file --------------------------------------------- 

TPM_values <- sleuth_to_matrix(so_rhox, 'obs_norm', 'tpm')
head(TPM_values)

TPM_values <- cbind(rownames(TPM_values),TPM_values)  #change rownames to a column
rownames(TPM_values) <- NULL  #remove rownames
colnames(TPM_values)[1] <- "target_id"  #rename the column as target_id,  later use this column to merge with
head(TPM_values)

merged_data <-merge(results_table,TPM_values,by=c("target_id"))
head(merged_data)

#--------------- Calculate the Fold Change from TPM-----------------------------
wildtype<- which(s2c$condition=="BV2_HA" )
overexpressed<- which(s2c$condition=="BV2_Rhox5")

control <- s2c[c(wildtype),]$sample
test <- s2c[c(overexpressed),]$sample

norm_test <-subset(so_rhox$obs_norm,so_rhox $obs_norm$sample %in% test)
norm_control <-subset(so_rhox $obs_norm,so_rhox $obs_norm$sample %in% control)
mean_t<-aggregate(tpm~target_id,data=norm_test,FUN=function(x) c(mean=mean(x)))
colnames(mean_t)<-c("target_id","BV2_Rhox_TPM")
head(mean_t)

mean_ctrl<-aggregate(tpm~target_id,data=norm_control,FUN=function(x) c(mean=mean(x)))
colnames(mean_ctrl)<-c("target_id","Control_TPM")
merged_mean <-merge(mean_t,mean_ctrl,by=c("target_id"))
merged_mean <- mutate(merged_mean , Ratio = (BV2_Rhox_TPM+1)/(Control_TPM+1), log2 = log2(Ratio))
tpm_results  <-left_join(merged_data ,merged_mean ,by=c("target_id"))
head(tpm_results )

#Remove redundant columns
drops <- c("ens_gene.y","ext_gene.y", "description.y")
write.table(tpm_results[ , !(names(tpm_results) %in% drops)] ,file=paste("sleuth_table_tx.tsv"), quote=FALSE, sep='\t', row.names = FALSE )

#-------------------------------------------------------------------------------

saveRDS(so_rhox, file = 'so.rds')
# Manually download Summaries -> Processed Data
sleuth_live(so_rhox)

### Plotting
## Maps
# PCA
plot_pca(so_rhox, color_by = 'condition', units = "est_counts")  
ggsave("pca_est.png", width = 190, height = 150, units ="mm")
plot_pca(so_rhox, color_by = 'condition', units = "tpm")  
ggsave("pca_tpm.png", width = 190, height = 150, units ="mm")
# Heatmap
png(filename = "heatmap.png",
    width = 190, height = 150, units = "mm", res=400)
plot_sample_heatmap(so_rhox)
dev.off()

## Analyses
# MA
plot_ma(so_rhox, test = "CHANGE_TO_YOUR_MODEL", test_type = "wt", which_model = "full",
        sig_level = 0.1)
ggsave("ma.png", width = 190, height = 150, units ="mm")
# Volcano
plot_volcano(so_rhox, test = "CHANGE_TO_YOUR_MODEL", test_type = "wt", which_model = "full",
             sig_level = 0.1, point_alpha = 0.2, sig_color = "red",
             highlight = NULL)
ggsave("volcano.png", width = 190, height = 150, units ="mm")

## Diagnostics
# Mean Variance Plot
plot_mean_var(so_rhox, which_model = "full", point_alpha = 0.4,
              point_size = 2, point_colors = c("black", "dodgerblue"),
              smooth_alpha = 1, smooth_size = 0.75, smooth_color = "red")
ggsave("mean_var.png", width = 190, height = 150, units ="mm")
# Quantile - Quantile
plot_qq(so_rhox, test = "CHANGE_TO_YOUR_MODEL", test_type = "wt", which_model = "full",
        sig_level = 0.1, point_alpha = 0.2, sig_color = "red",
        highlight = NULL, highlight_color = "green", line_color = "blue")
ggsave("qq.png", width = 190, height = 150, units ="mm")

## Other
#Plot technical variance versus observed variance
plot_vars(so_rhox, test = "CHANGE_TO_YOUR_MODEL", test_type = "wt", which_model = "full",
          sig_level = 0.1, point_alpha = 0.2, sig_color = "red", xy_line = TRUE,
          xy_line_color = "red", highlight = NULL, highlight_color = "green")

#Plot Loadings
plot_loadings(so_rhox, use_filtered = TRUE, sample = NULL, pc_input = NULL,
              units = "est_counts", pc_count = NULL, scale = FALSE,
              pca_loading_abs = TRUE)



