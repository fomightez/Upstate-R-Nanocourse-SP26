# install packages to R library

install.packages("data.table") # data loader
install.packages("plyr") # data management library
install.packages("BiocManager") # access to bioconductor library
install.packages("ggplot2") # plotting 
install.packages("reshape2") # changing shape of data frames
BiocManager::install("limma") # statistical tools for transcriptome data
BiocManager::install("edgeR") # processing & analysis of rna-seq counts
BiocManager::install("GEOquery") # download geo data directly into r environment
BiocManager::install("sva")
install.packages("cowplot") # combine plots into single figure
install.packages("RColorBrewer") # color palettes for figures

# load packages into environment
require(data.table)
require(plyr)
require(ggplot2)
require(reshape2)
require(limma)
require(edgeR)
require(GEOquery)
require(sva)
require(cowplot)
require(RColorBrewer)

# Download RNA-seq count data
gse = getGEOSuppFiles("GSE124326", fetch_files = TRUE) 
print(gse)

# Download sample meta-data
gse = getGEO(GEO = "GSE124326", destdir = "~/Downloads")
pheno <- pData(gse[[1]])  # import sample meta-data

# load RNA-seq counts
exprs = data.table::fread("/Users/hessjo/GSE124326/GSE124326_count_matrix.txt.gz")
dim(exprs) # 57820 rows, 481 columns

# set rna-seq as DGEList object from edgeR
edat = edgeR::DGEList(exprs)

# estimate library-size normalization factors (trimmed mean of M-values or TMM)
edat = edgeR::calcNormFactors(edat, method = "TMM")

# counts per million transformation
cpm = edgeR::cpm(edat, log = "FALSE")

# filter out lowly expressed transcripts
n_samples = ncol(cpm)

keep_transcripts = which(rowSums(cpm >= 1) >= (0.5 * n_samples))

edat_filter = edat[keep_transcripts, ,keep.lib.sizes = FALSE]

# log2 cpm transformation of retained transcripts
log2cpm <- edgeR::cpm(
  edat_filter,
  log = TRUE, # log2 transform
  prior.count = 1 # constant to add to values to avoid taking log of zero
)
rownames(log2cpm) = edat_filter$genes$gene

# inspect log2cpm distributions across 25 random samples
sample_index = sample(1:ncol(log2cpm), 25)
boxplot(log2cpm[,sample_index], las = 2)  # base R plotting method

# alternate: plotting log2-cpm distributions with ggplot2
long = reshape2::melt(log2cpm[,sample_index]) # first, convert wide format to long
long$Var2 = gsub("[.counts]", "", long$Var2) # long format of data

# then make plot with ggplot2
log2cpm_plot = ggplot(long, aes(x = Var2, y = value)) + # sets up basic plot aesthetics 
  geom_boxplot(fill = 'salmon1', outlier.color = 'grey') + # plot choice
  theme_bw() + # plot "theme" (black and white)
  xlab("Sample name") + # label x-axis
  ylab(expression(paste("log"[2], "(CPM)"))) + # label y-axis
  theme(panel.grid = element_blank(), # remove unnecessary gridlines from plot
        axis.text = element_text(size=11, color='black'), # increase text size for x- and y- axes
        axis.text.x = element_text(angle=45, hjust=1)) # tilt the x-axis labels on 45-degree angle

print(log2cpm_plot)
 
# inspect reduced dimensions
pca = prcomp(t(log2cpm)) # run principal component analysis 
dim(pca$x) # samples x samples eigenvectors

pca_df = data.frame(pca$x)

var_exp = (pca$sdev^2)/sum(pca$sdev^2) # proportion of variance explained by individual components
x_var_label = paste0("PC1 (",round(var_exp[[1]]*100,1), "%)")
y_var_label = paste0("PC2 (",round(var_exp[[2]]*100,1), "%)")
z_var_label = paste0("PC3 (",round(var_exp[[3]]*100,1), "%)")

rownames(pca_df) = gsub("[.counts]", "", rownames(pca_df))
identical(rownames(pca_df), pheno$title)

# pc1 vs pc2
viz_1 = ggplot(pca_df, aes(x = PC1, y = PC2, fill = pheno$characteristics_ch1.8)) + # set up basic plot aesthetics 
  geom_point(pch=21, stroke=0.2, size=2.5) + # make scatterplot 
  scale_fill_brewer("Sequencing plate", palette = 'Set1') +
  xlab(x_var_label) + # x axis label
  ylab(y_var_label) # y axis label

# pc1 vs pc3
viz_2 = ggplot(pca_df, aes(x = PC1, y = PC3, fill = pheno$characteristics_ch1.8)) + # set up basic plot aesthetics 
  geom_point(pch=21, stroke=0.2, size=2.5) + # make scatterplot 
  scale_fill_brewer("Sequencing plate", palette = 'Set1') +
  xlab(x_var_label) + # x axis label
  ylab(z_var_label) # y axis label

# pc2 vs pc3
viz_3 = ggplot(pca_df, aes(x = PC2, y = PC3, fill = pheno$characteristics_ch1.8)) + # set up basic plot aesthetics 
  geom_point(pch=21, stroke=0.2, size=2.5) + # make scatterplot 
  scale_fill_brewer("Sequencing plate", palette = 'Set1') +
  xlab(y_var_label) + # x axis label
  ylab(z_var_label) # y axis label


# batch correction looks justifiable (include primary factor + common sources of gene expression variance)
pheno$`age:ch1` = as.numeric(pheno$`age:ch1`)
pheno$`Sex:ch1` = factor(pheno$`Sex:ch1`)
pheno$`rin:ch1` = as.numeric(pheno$`rin:ch1`)
pheno$`lithium use (non-user=0, user = 1):ch1` = factor(pheno$`lithium use (non-user=0, user = 1):ch1`)
table(pheno$`bipolar disorder diagnosis:ch1`) # diagnostic groups
pheno$`bipolar disorder diagnosis:ch1` = factor(pheno$`bipolar disorder diagnosis:ch1`, levels=c("Control","BP1","BP2"))

# missing age info?
missing_age = which(is.na(pheno$`age:ch1`))
pheno_sub = pheno[-missing_age, ]
log2cpm_sub = log2cpm[,-missing_age]

pheno_sub_df = data.frame(age = pheno_sub$`age:ch1`,
                          sex = pheno_sub$`Sex:ch1`,
                          rin = pheno_sub$`rin:ch1`,
                          lithium = pheno_sub$`lithium use (non-user=0, user = 1):ch1`,
                          dx = pheno_sub$`bipolar disorder diagnosis:ch1`,
                          batch = pheno_sub$`sequencing plate:ch1`)

design = model.matrix(~age + sex + rin + dx, data = pheno_sub_df)

combat = sva::ComBat(dat = log2cpm_sub, batch = pheno_sub_df$batch, mod = design)


# repeat pca but with combat-corrected data
pca = prcomp(t(combat)) # run principal component analysis 
dim(pca$x) # samples x samples eigenvectors

pca_df = data.frame(pca$x)

var_exp = (pca$sdev^2)/sum(pca$sdev^2) # proportion of variance explained by individual components
x_var_label = paste0("PC1 (",round(var_exp[[1]]*100,1), "%)")
y_var_label = paste0("PC2 (",round(var_exp[[2]]*100,1), "%)")
z_var_label = paste0("PC3 (",round(var_exp[[3]]*100,1), "%)")
 
# manually start Set1 palette at 2nd (1st color [red] was NA in first set of plots, so ignore that color)
set1_colors = RColorBrewer::brewer.pal(n = 9, name = "Set1")[2:9]

# pc1 vs pc2
batch_viz_1 = ggplot(pca_df, aes(x = PC1, y = PC2, fill = pheno_sub$characteristics_ch1.8)) + # set up basic plot aesthetics 
  geom_point(pch=21, stroke=0.2, size=2.5) + # make scatterplot 
  scale_fill_manual("Sequencing plate", values = set1_colors) +
  xlab(x_var_label) + # x axis label
  ylab(y_var_label) # y axis label

# pc1 vs pc3
batch_viz_2 = ggplot(pca_df, aes(x = PC1, y = PC3, fill = pheno_sub$characteristics_ch1.8)) + # set up basic plot aesthetics 
  geom_point(pch=21, stroke=0.2, size=2.5) + # make scatterplot 
  scale_fill_manual("Sequencing plate", values = set1_colors) +
  xlab(x_var_label) + # x axis label
  ylab(z_var_label) # y axis label

# pc2 vs pc3
batch_viz_3 <- ggplot(pca_df, aes(x = PC2, y = PC3, fill = pheno_sub$characteristics_ch1.8)) +
  geom_point(pch=21, stroke=0.2, size=2.5) +
  scale_fill_manual("Sequencing plate", values = set1_colors) +
  xlab(y_var_label) +
  ylab(z_var_label)


# combine all batch-visualiations into single figure
pre_combat = cowplot::plot_grid(viz_1, viz_2, viz_3, ncol = 1)
post_combat = cowplot::plot_grid(batch_viz_1, batch_viz_2, batch_viz_3, ncol = 1)
cowplot::plot_grid(pre_combat, post_combat, ncol = 2)