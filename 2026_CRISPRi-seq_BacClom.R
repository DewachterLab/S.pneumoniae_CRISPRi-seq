# Script used for analysis of CRISPRi-seq results
# S. pneumoniae D39V prs1::lacI bgaA:Plac-dcas9 ZIP::P3-gRNA library
# Grown for 14 generations +/- IPTG +/- Bac 0.5 µg/ml +/- Clom 4 µg/ml
# 4 replicates in each condition
# Illumina sequencing > 2FAST2Q > DESeq2 (= current script)


# Preliminaries -----------------------------------------------------------

library(tidyverse)
library(readxl)
library(DESeq2)
library(ggbeeswarm)
library(data.table)
library(apeglm)
library(ggtext)
library(plotly)
library(ggrepel)
library(writexl)
library(ggtext)

# Choose correct working directory
setwd("")


# Data import & wrangling -------------------------------------------------

# Load read counts from 2FAST2Q compiled output
raw <- read.csv2("compiled_Processed.csv")
# Turn first column of read counts (= gRNA numbers) into row names
cts <- as.matrix(raw[,-1]) # select counts only
rownames(cts) <- raw[,1] # add gRNA row names


# Load column data file which has sample information
coldata <- read_excel("coldata.xlsx")
# Add sample number as row names for each condition in column data file
rownames(coldata) <- colnames(cts)




# Associating gRNA with target gene names ---------------------------------
targetdf <- read.csv(file = "target_list_ldw.csv")



# Create DESeq2 object ----------------------------------------------------

# Create DESeq2 object with counts and sample info
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~  Condition * Induction)





# Plotting read counts (raw & normalized) ---------------------------------

# Put counts table in tidy format (every row is one data point)
cts_tidy <- raw %>% pivot_longer(cols = !c(X.sgRNA), names_to = "Sample", values_to = "Counts")
cts_tidy <- data.frame(cts_tidy,
                       coldata[match(cts_tidy$Sample, 
                                     rownames(coldata)), ])

# PLOT OF RAW READ COUNTS
ggplot(cts_tidy, aes(x = Condition, y = Counts, col = Induction, group = Sample)) + 
  geom_quasirandom(dodge.width = 1) + 
  scale_color_brewer(palette = "Dark2") + # colorblind-friendly palette
  labs(x = "Condition", 
       y = "Raw counts", 
       col = "CRISPRi induction") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        axis.ticks.length = unit(0.2,"cm"),
        axis.text.x = element_text(size=16,color="black"), 
        axis.text.y = element_text(size=16,color="black"),
        axis.title.x = element_text(size=18,margin=margin(20,0,0,0),color="black"),
        axis.title.y = element_text(size=18,margin=margin(0,20,0,0),color="black"),
        legend.title = element_text(size=18,color="black"),
        legend.text = element_text(size=16,color="black"),
        plot.title = element_text(size=18, color="black",hjust=0.5))


# Normalize read counts using DESeq2's method
dds <- estimateSizeFactors(dds)
cts_norm <- counts(dds, normalized = TRUE)

# Turn cts_norm into a df which is needed for pivot_longer
cts_norm_df <- bind_cols(rownames(cts_norm), as.data.frame(cts_norm)) 
names(cts_norm_df)[names(cts_norm_df) == "...1"] <- "X.sgRNA"

# Put normalized counts table in tidy format (every row is one data point)
cts_norm_tidy <- cts_norm_df %>% pivot_longer(cols = !c(X.sgRNA), names_to = "Sample", values_to = "Counts")
cts_norm_tidy <- data.frame(cts_norm_tidy,
                       coldata[match(cts_norm_tidy$Sample, 
                                     rownames(coldata)), ])

# PLOT OF NORMALIZED READ COUNTS
ggplot(cts_norm_tidy, aes(x = Condition, y = Counts, col = Induction, group = Sample)) + 
  geom_quasirandom(dodge.width = 1) + 
  scale_color_brewer(palette = "Dark2") + # colorblind-friendly palette
  labs(x = "Condition", 
       y = "Raw counts", 
       col = "CRISPRi induction") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        axis.ticks.length = unit(0.2,"cm"),
        axis.text.x = element_text(size=16,color="black"), 
        axis.text.y = element_text(size=16,color="black"),
        axis.title.x = element_text(size=18,margin=margin(20,0,0,0),color="black"),
        axis.title.y = element_text(size=18,margin=margin(0,20,0,0),color="black"),
        legend.title = element_text(size=18,color="black"),
        legend.text = element_text(size=16,color="black"),
        plot.title = element_text(size=18, color="black",hjust=0.5))



# Principal component analysis (PCA) --------------------------------------

rld <- rlog(dds, blind = FALSE)
pcaData <- plotPCA(rld, intgroup = c("Condition", "Induction"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = Condition, shape = Induction)) +
  geom_point(size = 10) +
  scale_color_brewer(palette = "Dark2") + 
  scale_shape_manual(values = c("-", "+")) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  guides(col = guide_legend(override.aes = list(size = 4))) + 
  coord_fixed() +
  theme_bw(base_size = 18) +
  theme(axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        axis.ticks.length = unit(0.2,"cm"),
        axis.text.x = element_text(size=16,color="black"), 
        axis.text.y = element_text(size=16,color="black"),
        axis.title.x = element_text(size=18,margin=margin(20,0,0,0),color="black"),
        axis.title.y = element_text(size=18,margin=margin(0,20,0,0),color="black"),
        legend.title = element_text(size=18,color="black"),
        legend.text = element_text(size=16,color="black"),
        plot.title = element_text(size=18, color="black",hjust=0.5))




# Differential fitness analysis -------------------------------------------
alpha <- 0.05 # Significance threshold
log2FCtest <- 1 # Threshold for meaningful depletion/enrichment (1 = 2-fold change)

dds <- DESeq(dds)
resultsNames(dds)


# Comparing + vs - IPTG in condition without AB
# Extract results
minAB_res <- results(dds, name = "Induction_plusIPTG_vs_minIPTG",
                     lfcThreshold = log2FCtest, alpha = alpha)
# View summary of results
summary(minAB_res)
# shrink L2FC to make more robust
minAB_res_shr <- lfcShrink(dds, res = minAB_res, type = "apeglm",
                           coef = "Induction_plusIPTG_vs_minIPTG")
# add targets of gRNAs
minAB_res_shr$target <- targetdf$target_names





# Comparing + vs - IPTG in condition with Bac
# Set plusBac as reference
dds_plusBac_ref <- dds
dds_plusBac_ref$Condition <- relevel(dds_plusBac_ref$Condition,
                                     ref = "plusBac")
resultsNames(dds_plusBac_ref)
# rerun the hypothesis testing part of DESeq()
dds_plusBac_ref <- nbinomWaldTest(dds_plusBac_ref)
# extracting results exactly like before
plusBac_res <- results(dds_plusBac_ref, name = "Induction_plusIPTG_vs_minIPTG", 
                       lfcThreshold = log2FCtest, alpha = alpha)
summary(plusBac_res)
# shrink L2FC to make more robust
plusBac_res_shr <- lfcShrink(dds_plusBac_ref, res = plusBac_res, type = "apeglm",
                             coef = "Induction_plusIPTG_vs_minIPTG")
# add targets of gRNAs
plusBac_res_shr$target <- targetdf$target_names





# Comparing + vs - IPTG in condition with Clom
# Set plusClom as reference
dds_plusClom_ref <- dds
dds_plusClom_ref$Condition <- relevel(dds_plusClom_ref$Condition,
                                      ref = "plusClom")
resultsNames(dds_plusClom_ref)
# rerun the hypothesis testing part of DESeq()
dds_plusClom_ref <- nbinomWaldTest(dds_plusClom_ref)
# extracting results exactly like before
plusClom_res <- results(dds_plusClom_ref, name = "Induction_plusIPTG_vs_minIPTG", 
                        lfcThreshold = log2FCtest, alpha = alpha)
summary(plusClom_res)
# shrink L2FC to make more robust
plusClom_res_shr <- lfcShrink(dds_plusClom_ref, res = plusClom_res, type = "apeglm",
                              coef = "Induction_plusIPTG_vs_minIPTG")
# add targets of gRNAs
plusClom_res_shr$target <- targetdf$target_names





# INTERACTION: Comparing L2FC caused by IPTG induction for + vs - Bacitracin
# Extract results of the interaction term
int_res_Bac <- results(dds, name = "ConditionplusBac.InductionplusIPTG",
                       lfcThreshold = log2FCtest, alpha = alpha)

# View summary of results
summary(int_res_Bac)
# shrink L2FC to make more robust
int_res_Bac_shr <- lfcShrink(dds, res = int_res_Bac, type = "apeglm",
                             coef = "ConditionplusBac.InductionplusIPTG")
# add targets of gRNAs
int_res_Bac_shr$target <- targetdf$target_names


# Make correlation plot
# Point size scales with adjusted P value
plusBac_intscat <- data.frame("Tag" = rownames(plusBac_res_shr), 
                              "Target" = targetdf$target_names, 
                              "log2FCplusBac" = plusBac_res_shr$log2FoldChange, 
                              "SElog2FCplusBac" = plusBac_res_shr$lfcSE, 
                              "log2FCminAB" = minAB_res_shr$log2FoldChange, 
                              "SElog2FCminAB" = minAB_res_shr$lfcSE, 
                              "log2FCint" = int_res_Bac_shr$log2FoldChange, 
                              "SElog2FCint" = int_res_Bac_shr$lfcSE, 
                              "padj_int" = int_res_Bac_shr$padj, 
                              "Interaction" = ifelse(int_res_Bac_shr$padj < alpha, 
                                                     "yes", "no"))

gg_intscat_Bac_size <- ggplot(plusBac_intscat, 
                              aes(x = log2FCminAB, 
                                  y = log2FCplusBac, 
                                  label = ifelse(abs(log2FCint) > 2^log2FCtest & Interaction == "yes", 
                                                 Target, NA), 
                                  size = -log10(padj_int), 
                                  col = Interaction)) + 
  geom_vline(xintercept = c(0), lty = 2, col = "grey", linewidth = 1) +
  geom_hline(yintercept = c(0), lty = 2, col = "grey", linewidth = 1) + 
  geom_abline(intercept = c(0), slope = 1, lty = 2, col = "grey", linewidth = 1) + 
  geom_point() + 
  scale_color_brewer(palette = "Dark2", na.translate = FALSE) +
  scale_size_continuous(range = c(1, 6), 
                        breaks = c(0, 20, 40, 60),
                        limits = c(0, 60)) +
  xlim(-14, 2) +
  ylim(-14, 6) +
  labs(x = paste("Fitness without Bac (Log<sub>2</sub>FC)"), 
       y = paste("Fitness with Bac (Log<sub>2</sub>FC)"), 
       size = expression(-Log["10"](P["adj"]))) + 
  theme_classic(base_size = 18) + 
  theme(axis.title.x = element_markdown(), 
        axis.title.y = element_markdown(),
        #legend.position = "none",
        axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        axis.ticks.length = unit(0.2,"cm"),
        axis.text.x = element_text(size=16,color="black"), 
        axis.text.y = element_text(size=16,color="black"),
        legend.title = element_text(size=18,color="black"),
        legend.text = element_text(size=16,color="black"))
gg_intscat_Bac_size
ggplotly(gg_intscat_Bac_size)




# INTERACTION: Comparing L2FC caused by IPTG induction for + vs - Clomiphene
# Extract results of the interaction term
int_res_Clom <- results(dds, name = "ConditionplusClom.InductionplusIPTG",
                        lfcThreshold = log2FCtest, alpha = alpha)
# View summary of results
summary(int_res_Clom)
# shrink L2FC to make more robust
int_res_Clom_shr <- lfcShrink(dds, res = int_res_Clom, type = "apeglm",
                              coef = "ConditionplusClom.InductionplusIPTG")
# add targets of gRNAs
int_res_Clom_shr$target <- targetdf$target_names



# make correlation plot
# Point size scales with adjusted P value
plusClom_intscat <- data.frame("Tag" = rownames(plusClom_res_shr), 
                              "Target" = targetdf$target_names, 
                              "log2FCplusCom" = plusClom_res_shr$log2FoldChange, 
                              "SElog2FCplusCom" = plusClom_res_shr$lfcSE, 
                              "log2FCminAB" = minAB_res_shr$log2FoldChange, 
                              "SElog2FCminAB" = minAB_res_shr$lfcSE, 
                              "log2FCint" = int_res_Clom_shr$log2FoldChange, 
                              "SElog2FCint" = int_res_Clom_shr$lfcSE, 
                              "padj_int" = int_res_Clom_shr$padj, 
                              "Interaction" = ifelse(int_res_Clom_shr$padj < alpha, 
                                                     "yes", "no"))

gg_intscat_Clom_size <- ggplot(plusClom_intscat, 
                               aes(x = log2FCminAB, 
                                   y = log2FCplusCom, 
                                   label = ifelse(abs(log2FCint) > 2^log2FCtest & Interaction == "yes", 
                                                  Target, NA), 
                                   size = -log10(padj_int), 
                                   col = Interaction)) + 
  geom_vline(xintercept = c(0), lty = 2, col = "grey", linewidth = 1) +
  geom_hline(yintercept = c(0), lty = 2, col = "grey", linewidth = 1) + 
  geom_abline(intercept = c(0), slope = 1, lty = 2, col = "grey", linewidth = 1) + 
  geom_point() + 
  scale_color_brewer(palette = "Dark2", na.translate = FALSE) + 
  scale_size_continuous(range = c(1, 6), 
                        breaks = c(0, 20, 40, 60),
                        limits = c(0, 60)) +
  xlim(-14, 2) +
  ylim(-14, 6) +
  labs(x = paste("Fitness without Clom (Log<sub>2</sub>FC)"), 
       y = paste("Fitness with Clom (Log<sub>2</sub>FC)"), 
       size = expression(-Log["10"](P["adj"]))) + 
  theme_classic(base_size = 18) + 
  theme(axis.title.x = element_markdown(), 
        axis.title.y = element_markdown(),
        #legend.position = "none",
        axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        axis.ticks.length = unit(0.2,"cm"),
        axis.text.x = element_text(size=16,color="black"), 
        axis.text.y = element_text(size=16,color="black"),
        legend.title = element_text(size=18,color="black"),
        legend.text = element_text(size=16,color="black"))
gg_intscat_Clom_size
ggplotly(gg_intscat_Clom_size)



# INTERACTION: Comparing L2FC caused by IPTG induction for Bacitracin vs Clomiphene
# Extract results of the interaction term
int_res_BvC <- results(dds_plusClom_ref, name = "ConditionplusBac.InductionplusIPTG",
                       lfcThreshold = log2FCtest, alpha = alpha)
# View summary of results
summary(int_res_BvC)
# shrink L2FC to make more robust
int_res_BvC_shr <- lfcShrink(dds, res = int_res_BvC, type = "apeglm",
                             coef = "ConditionplusBac.InductionplusIPTG")
# add targets of gRNAs
int_res_BvC_shr$target <- targetdf$target_names


# make interaction plat
# Point size scales with adjusted P value
plusBvC_intscat <- data.frame("Tag" = rownames(plusClom_res_shr), 
                              "Target" = targetdf$target_names, 
                              "log2FCplusBac" = plusBac_res_shr$log2FoldChange, 
                              "SElog2FCplusBac" = plusBac_res_shr$lfcSE, 
                              "log2FCplusClom" = plusClom_res_shr$log2FoldChange, 
                              "SElog2FCplusClom" = plusClom_res_shr$lfcSE, 
                              "log2FCint" = int_res_BvC_shr$log2FoldChange, 
                              "SElog2FCint" = int_res_BvC_shr$lfcSE, 
                              "padj_int" = int_res_BvC_shr$padj, 
                              "Interaction" = ifelse(int_res_BvC_shr$padj < alpha, 
                                                     "yes", "no"))

gg_intscat_BvC_size <- ggplot(plusBvC_intscat, 
                              aes(x = log2FCplusClom, 
                                  y = log2FCplusBac, 
                                  label = ifelse(abs(log2FCint) > 2^log2FCtest & Interaction == "yes", 
                                                 Target, NA), 
                                  size = -log10(padj_int), 
                                  col = Interaction)) + 
  geom_vline(xintercept = c(0), lty = 2, col = "grey", linewidth = 1) +
  geom_hline(yintercept = c(0), lty = 2, col = "grey", linewidth = 1) + 
  geom_abline(intercept = c(0), slope = 1, lty = 2, col = "grey", linewidth = 1) + 
  geom_point() + 
  scale_color_brewer(palette = "Dark2", na.translate = FALSE) + 
  scale_size_continuous(range = c(1, 6), 
                        breaks = c(0, 20, 40, 60),
                        limits = c(0, 100)) +
  xlim(-14, 2) +
  ylim(-14, 6) +
  labs(x = paste("Fitness with Clom (Log<sub>2</sub>FC)"), 
       y = paste("Fitness with Bac (Log<sub>2</sub>FC)"), 
       size = expression(-Log["10"](P["adj"]))) + 
  theme_classic(base_size = 18) + 
  theme(axis.title.x = element_markdown(), 
        axis.title.y = element_markdown(),
        #legend.position = "none",
        axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        axis.ticks.length = unit(0.2,"cm"),
        axis.text.x = element_text(size=16,color="black"), 
        axis.text.y = element_text(size=16,color="black"),
        legend.title = element_text(size=18,color="black"),
        legend.text = element_text(size=16,color="black"))
gg_intscat_BvC_size
ggplotly(gg_intscat_BvC_size)

