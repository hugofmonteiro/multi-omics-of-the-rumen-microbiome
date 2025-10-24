# Diff gene expression and volcano plots
library(DESeq2)

setwd("~/Downloads/")

# Adjust the path to your files; countData was also performed with genome-guided and microbial isoform counts
countData <- read.csv("mRNA_only_byIsoform_counts_matrix.csv", sep = "\t", header = TRUE, row.names = 1) 
colData <- read.csv("metadata.csv", row.names=1)

colnames(countData) <- gsub("X", "", colnames(countData)) 
colnames(countData) <- gsub("\\.", "-", colnames(countData)) 
colnames(countData) <- gsub("R-", "XR-", colnames(countData)) 

#finding mismatches
mismatched <- setdiff(rownames(colData), colnames(countData))
print(mismatched)

# Check again to ensure they match
head(colnames(countData))
head(rownames(colData))

# Reorder the columns of countData to match the order of rownames in colData
countData <- countData[, rownames(colData)]
countData <- round(countData)
countData <- as.matrix(countData)
storage.mode(countData) <- "integer"

# minimum count for a transcript to be considered expressed
min_count <- 50
# minimum number of samples in which a transcript needs to be expressed (decision made based on the diversity of the microbiome and sparsity of transcript abundances)
min_samples <- 1

# Splitting count data based on 'Group' levels
countData_least <- countData[, colData$Group == 'Least']
countData_most <- countData[, colData$Group == 'Most']

# Applying filtering criteria
keep_transcripts_least <- rowSums(countData_least >= min_count) >= min_samples
keep_transcripts_most <- rowSums(countData_most >= min_count) >= min_samples

# Combine results: keep a gene if it meets the criteria in either subset
keep_transcripts <- keep_transcripts_least | keep_transcripts_most

# Filter the original count data
countData_filtered <- countData[keep_transcripts, ]

# DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = countData_filtered,
                              colData = colData,
                              design = ~ Group)
dds <- DESeq(dds)
res <- results(dds)
summary(res, alpha=0.05)

# converting the results to data.frame for plotting
results <- as.data.frame(res)

# saving datasets for sanity check
write.csv(as.data.frame(results), file="DESeq2_results.csv") # also for genome-guided and microbial transcripts

library(ggplot2)

# Define the significance threshold
pvalue_threshold <- 0.05
log2FC_threshold <- 1

significant_results <- subset(results, pvalue < pvalue_threshold & abs(log2FoldChange) > log2FC_threshold)
write.csv(as.data.frame(significant_results), file="DESeq2_significantTranscripts.csv")

# I had to open DESeq2_significantTranscripts.csv and remove the transcripts with -15 < L2FC > 15 so most transcripts could be displayed
# Uploading cleaned dataset without extreme L2FC to join plots later
significant_results_LFCbelow15 <- read.csv("DESeq2_significantTranscripts_LFCbelow15.csv", row.names=1)
# Once I had a volcano with all transcripts and another without extreme L2FC, I manually joined both plots in InkScape.

# Adding labels to significantly different with FDR transcripts
significant_results_LFCbelow15$status <- ifelse(significant_results_LFCbelow15$padj < pvalue_threshold & significant_results_LFCbelow15$log2FoldChange >  0, "Up",
                                                ifelse(significant_results_LFCbelow15$padj < pvalue_threshold & significant_results_LFCbelow15$log2FoldChange <  0, "Down", "NS"))
significant_results_LFCbelow15$status <- factor(significant_results_LFCbelow15$status, levels = c("NS","Down","Up"))

####################################
significant_results_LFCbelow15$logp <- -log10(significant_results_LFCbelow15$padj)
significant_results_LFCbelow15 <- significant_results_LFCbelow15[is.finite(significant_results_LFCbelow15$logp), ]

# Debugging the range of the data for both plots
cat("Range of log2FoldChange: ", range(significant_results_LFCbelow15$log2FoldChange, na.rm = TRUE), "\n")
cat("Range of logp: ", range(significant_results_LFCbelow15$logp, na.rm = TRUE), "\n")

# Normalize baseMean for point size scaling
significant_results_LFCbelow15$size_scaled <- ifelse(significant_results_LFCbelow15$baseMean <= 100, 2,
                                                     ifelse(significant_results_LFCbelow15$baseMean <= 500, 4,
                                                            ifelse(significant_results_LFCbelow15$baseMean <= 1000, 6,
                                                                   10)))

# volcano plot using ggplot
volcano_plot <- ggplot(significant_results_LFCbelow15, aes(x = log2FoldChange, y = logp, size = size_scaled, fill= status)) +
  geom_point(shape = 21, color = "darkgrey", alpha = 0.5) +
  scale_size(range = c(4, 10), name = "BaseMean (scaled)") +
  scale_fill_manual(values = c(NS = "#DDDDDD", Down = "#234da0", Up = "#ff0000")) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )+
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  xlab(expression(Log[2]*" Fold Change")) +
  ylab(expression(-log[10]*" (Adj P-value)"))
volcano_plot

ggsave("Volcano_Plot_DiffExpression_LFCbelow15.svg", plot = volcano_plot, dpi = 600)
