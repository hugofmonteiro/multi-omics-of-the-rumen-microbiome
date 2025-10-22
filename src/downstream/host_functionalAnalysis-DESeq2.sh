#Remove Stored Variables
rm(list = ls())

'%!in%' <- function(x,y)!('%in%'(x,y))

# if we're using Rstudio, set the working directory to wherever this file is saved
if(.Platform$GUI == "RStudio"){
  if( TRUE %!in% grepl(pattern = "rstudioapi", x = rownames(installed.packages())) ){
    install.packages('rstudioapi')
  }
  library(rstudioapi)
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
} else {
	setwd("~/storage/ctbrownroot/Hugo/BOSTAURUS/HISAT2/ALIGN/")
}

# Package names
packages <- c("viridis", "ggVennDiagram", "ggrepel","plyr","ggpubr","DESeq2","purrr", "tidyverse", "tidyr","dplyr","stringr","ggplot2", "reshape2", "viridis", "ggforce", "BiocManager", "openxlsx","DGEobj.utils", "edgeR")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

list.files(path = "/home/jknorris/storage/ctbrownroot/Hugo/BOSTAURUS/HISAT2/ALIGN", pattern = "*count")

Metadata <- read.csv( "metadata_hugo.csv")
# set up a table describing your samples, variables that you want to compare, etc.
#directory <- paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
directory <- "/home/jknorris/storage/ctbrownroot/Hugo/BOSTAURUS/HISAT2/ALIGN/"
sampleFiles <- grep("*.counts$",list.files(directory),value=TRUE)

#list.files(path = ".", pattern = "*.count")

##############################
Metadata$files <- sampleFiles

sampleTable <- data.frame(sampleName = Metadata$X.IDSEQ,
                          fileName = Metadata$files,
                          condition = Metadata$RFI)

sampleTable$condition <- factor(sampleTable$condition)


  tmp_padj = 0.05

  directory <- getwd()

  rownames(Metadata) <- Metadata[,3]

  Metadata$RFI <- factor(Metadata$RFI)
  full_counts <- grep("*counts$",list.files(path=directory, full.names=TRUE),value=TRUE)
  full_counts <- paste(getwd(),Metadata$files, sep="/")

  read_in_feature_counts <- function(file){
    cnt <- read_tsv(file, col_names =T, comment = "#")
    cnt <- cnt %>% dplyr::select(-Chr, -Start, -End, -Strand, -Length)
    return(cnt)
  }

  raw_counts <- map(full_counts, read_in_feature_counts)
  raw_counts_df <- purrr::reduce(raw_counts, inner_join)
  raw_counts_df <- as.data.frame(raw_counts_df)
  rownames(raw_counts_df) <- raw_counts_df[,1]
  raw_counts_df <- raw_counts_df[,-1]
  rownames(Metadata) <- gsub(patter="18048XR-81-|.counts", replacement="", x=rownames(Metadata))

  colnames(Metadata) <- c("ID","GROUP","SampleFiles")
  colnames(raw_counts_df) <- gsub("/home/jknorris/storage/ctbrownroot/Hugo/BOSTAURUS//HISAT2/ALIGN/|.sam","", colnames(raw_counts_df))
  write.csv(x = raw_counts_df, file = "250113_Btaurus_mostvleast_rawcounts.csv", row.names = TRUE, col.names = TRUE)
  Metadata$GROUP <- factor(Metadata$GROUP)
rownames(Metadata) <- NULL
  deseq_featurecounts_Group <- DESeqDataSetFromMatrix(countData = raw_counts_df,
                                                      colData = Metadata,
                                                      design = ~ GROUP)

  deseq_featurecounts_Group <- estimateSizeFactors(deseq_featurecounts_Group)
deseq_featurecounts_Group <- estimateDispersionsGeneEst(deseq_featurecounts_Group)
  dispersions(deseq_featurecounts_Group) <- mcols(deseq_featurecounts_Group)$dispGeneEst

  DESeqResults_featurecounts_Group <- DESeq(deseq_featurecounts_Group)
  resultsNames(DESeqResults_featurecounts_Group)
  res_DESeqResults_Group_MvL <- as.data.frame(results(DESeqResults_featurecounts_Group, name="GROUP_Most_vs_Least", alpha=tmp_padj))

  openxlsx::write.xlsx(x = res_DESeqResults_Group_MvL, file = "DESEQ2_Group_Most_V_Least_HUGO.xlsx", rowNames=TRUE)

##################################################################################

res_deseq2 <- as.data.frame(res_DESeqResults_Group_MvL)

# Annotating DEGs based on their expression changes
res_deseq2 %>%
  mutate(Group = case_when(
    log2FoldChange >= 1 & padj <= 0.05 ~ "Up-regulated",
    log2FoldChange <= -1 & padj  <= 0.05 ~ "Down-regulated",
    TRUE ~ "Not-significant"
  )) -> res2

# Summarizing the count of DEGs in each category
table(res2$Group)

res2$logp <- -log10(res2$padj)


# Adding labels for significant genes in the volcano plot
res2$label <- ""
res2 <- res2[order(res2$padj),]
up.genes <- head(rownames(res2[which( (res2$log2FoldChange > 1) & (res2$padj < 0.05) ),]), 10)

down.genes <- head(rownames(res2[which((res2$log2FoldChange < -1) & (res2$padj < 0.05)),]), 10)
res2$label[match(c(up.genes, down.genes), rownames(res2))] <- c(up.genes, down.genes)

# Removing any rows with missing data
res2 <- na.omit(res2)

scatterplot <- ggscatter(res2, x="log2FoldChange",y="logp",
                              palette = c("#DDDDDD","#ff0000", "#234da0"),
                              size=3,
                              shape=21,
                              fill="Group",
                              color="grey",
                              alpha=0.5,
                              repel=T,
                              xlim=c(min(res2$log2FoldChange),max(res2$log2FoldChange)))+ theme_minimal() +
  geom_hline(yintercept=1.30,linetype="dashed")+
  geom_vline(xintercept= c(-1,1),linetype="dashed") +
  xlab(expression(Log[2]*" Fold Change")) +
  ylab(expression(-log[10]*" (Adj. P-value)"))

svg(filename = "Scatterplot.svg", width = 11, height = 11)
scatterplot
dev.off()

save.image(file="deseq.Rdata", compress=TRUE)
