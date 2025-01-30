# My Script 

#activate required libraries
library(tidyverse)
library(DESeq2)
library(readxl)
library(ggplot2)
library(ggrepel)
library(ggpubr)

#read in the gene counts and the sample information files

coldata <- read_excel("coldata.xlsx")
coldata$Participant <- as.factor(coldata$Participant)
gene_count <- read_excel("gene_count.xls")

#create a dds from the gene count and the sample information
dds_mRNA <- DESeqDataSetFromMatrix(
  countData = gene_count,
  colData = coldata,
  design = ~ Condition)

# convert the count data into a data frame.
gene_count <- as.data.frame(gene_count)


# need first to remove dublicate from the gene name columns, as it has to be unique.
gene_count2 <- gene_count [!duplicated(gene_count$gene_name),]

# convert the gene names column to row names and deltete the column woth gene name
gene_count3 <- gene_count2[,-1]
 # take gene name from gene_count2 and add it into gene_coun 3 in the row name
rownames(gene_count3) <- gene_count2 [,1]

# finally we can run the dds for gene_count3
dds_mRNA <- DESeqDataSetFromMatrix(
  countData = gene_count3,
  colData = coldata,
  design = ~ Condition)
# Filtration of count gene with very low counts.

low_expressed <- counts(dds_mRNA) >10
keep <- rowSums(low_expressed) == ncol(dds_mRNA)

Keep_dds_mRNA <- rowSums(counts(dds_mRNA)) >= 10
dds_mRNA <- dds_mRNA [keep,]

# PCA for gene counts
 # first we need to creeate a transform tool like the Log to reduce the variants
vsd <- vst (dds_mRNA)
plotPCA(vsd,intgroup = "Participant",ntop=500)
dds_mRNA_no_LPS <- dds_mRNA[,dds_mRNA$Condition!='LPS']
# how to see things in your file
unique(dds_mRNA$Condition)
unique(dds_mRNA_no_LPS$Condition)


# Repeat PCA without LPS
vsd2 <- vst (dds_mRNA_no_LPS)
plotPCA(vsd2,intgroup = 'Participant',ntop=500) 

# Plot PCA:
plotPCA (vsd2,intgroup = ("Participant"),ntop=500)

pcaData <- plotPCA(vsd2, intgroup=c("Condition", "Participant"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color= Participant , shape= Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

# to run the differential analysis:
#model comparing condition without adjustment
de_mRNA <- DESeq(dds_mRNA)
res_ctrl_vs_nLDL <- results(de_mRNA, 
       contrast = c("Condition", "CTRL","nLDL_NMVs"),tidy=TRUE) %>% arrange(pvalue)

# to see the rsults file:

head (res_ctrl_vs_nLDL)

# call result:

res_ctrl_vs_oxLDL<- results(de_mRNA, contrast = c("Condition", "CTRL", "oxLD_NMVs"),tidy=TRUE) %>%
  arrange(pvalue)

res_ctrl_vs_oxLDL

res_nLDL_vs_oxLDL<- results(de_mRNA, contrast = c("Condition", "nLDL_NMVs", "oxLD_NMVs"),tidy=TRUE) %>%
  arrange(pvalue)

res_nLDL_vs_oxLDL

resultsNames(de_mRNA)

# to export table:

write.csv(res_ctrl_vs_nLDL,'DeFF_CTRL vs nLDL_NMVs.csv')
write.csv(res_ctrl_vs_oxLDL,'DeFF_CTRL vs oxLDL_NMVs.csv')
write.csv(res_nLDL_vs_oxLDL,'DeFF_nLDL_NMVs vs nLDL_NMVs.csv')

#model comparing condition with adjustment for participant
dds_mRNA2 <- dds_mRNA
design(dds_mRNA2) <- ~ Condition + Participant

# to see the new design of the new file (2)
design(dds_mRNA2)


#model comparing condition with adjustment for Gender
dds_mRNA3 <- dds_mRNA
design(dds_mRNA3) <- ~ Condition + Participant

# Convert character variables to factors
# For example, if "condition" is a character variable
dds_mRNA$ Gender<- factor(dds_mRNA$Gender)

# to rerun the deseq after changing the design:
de_mRNA3 <- DESeq(dds_mRNA3)
design(de_mRNA3)

# call results:
res_ctrl_vs_nLDL3 <- results(de_mRNA3, contrast = c("Condition", "CTRL", "nLDL_NMVs"),tidy=TRUE) %>%
  arrange(pvalue)

# to see the results file:
head(res_ctrl_vs_nLDL3)

res_ctrl_vs_nLDL3 

# to create file for sig. only:
sig_res_ctrl_vs_nLDL3 <- subset(res_ctrl_vs_nLDL3, pvalue<0.05)

write.csv(sig_res_ctrl_vs_nLDL3, 'sig_diff_Pvalue_CTRL vs nLDL_NMVs3.csv')

# call result CTRL vs oxLDL:
res_ctrl_vs_oxLDL3<- results(de_mRNA3, contrast = c("Condition", "CTRL", "oxLD_NMVs"),tidy=TRUE) %>%
  arrange(pvalue)

# we can arrange by the P value when we call the results:
res_ctrl_vs_oxLDL3_by_Pvalue<- results(de_mRNA3, contrast = c("Condition", "CTRL", "oxLD_NMVs"),tidy=TRUE) %>%
  arrange(pvalue)

head(res_ctrl_vs_oxLDL3,20)

# to create file for sig. only:
sig_res_ctrl_vs_oxLDL3 <- subset(res_ctrl_vs_oxLDL3, pvalue<0.05)

write.csv(sig_res_ctrl_vs_oxLDL3, 'sig_diff_Pvalue CTRL vs oxLDL_NMVs3.csv')

# call result nLDL vs oxLDL:
res_nLDL_vs_oxLDL3<- results(de_mRNA3, contrast = c("Condition", "nLDL_NMVs", "oxLD_NMVs"),tidy=TRUE) %>%
  arrange(pvalue)

head(res_nLDL_vs_oxLDL3)


# to create file for sig. only:
sig_res_nLDL_vs_oxLDL3 <- subset(res_nLDL_vs_oxLDL3, pvalue<0.05)

write.csv(sig_res_nLDL_vs_oxLDL3, 'sig_diff_Pvalue nLDL-NMVs vs oxLDL_NMVs3.csv')

#volcano

#create a new coloumn called diffexpressed to categorise the expression to up, down and no

#ctrl vs nLDL
res_ctrl_vs_nLDL3$diffexpressed <- "NO"
res_ctrl_vs_nLDL3$diffexpressed[res_ctrl_vs_nLDL3$log2FoldChange > 0.0 & res_ctrl_vs_nLDL3$pvalue < 0.05] <- "UP"
res_ctrl_vs_nLDL3$diffexpressed[res_ctrl_vs_nLDL3$log2FoldChange < -0.0 & res_ctrl_vs_nLDL3$pvalue < 0.05] <- "DOWN"
head(res_ctrl_vs_nLDL3)
unique(res_ctrl_vs_nLDL3$diffexpressed)
table(res_ctrl_vs_nLDL3$diffexpressed) #up=349

#ctrl vs oxLDL
res_ctrl_vs_oxLDL3$diffexpressed <- "NO"
res_ctrl_vs_oxLDL3$diffexpressed[res_ctrl_vs_oxLDL3$log2FoldChange > 0.0 & res_ctrl_vs_oxLDL3$pvalue < 0.05] <- "UP"
res_ctrl_vs_oxLDL3$diffexpressed[res_ctrl_vs_oxLDL3$log2FoldChange < -0.0 & res_ctrl_vs_oxLDL3$pvalue < 0.05] <- "DOWN"
head(res_ctrl_vs_oxLDL3)
unique(res_ctrl_vs_oxLDL3$diffexpressed)
table(res_ctrl_vs_oxLDL3$diffexpressed)

#nLDL vs oxLDL
res_nLDL_vs_oxLDL3$diffexpressed <- "NO"
res_nLDL_vs_oxLDL3$diffexpressed[res_nLDL_vs_oxLDL3$log2FoldChange > 0.0 & res_nLDL_vs_oxLDL3$pvalue < 0.05] <- "UP"
res_nLDL_vs_oxLDL3$diffexpressed[res_nLDL_vs_oxLDL3$log2FoldChange < -0.0 & res_nLDL_vs_oxLDL3$pvalue < 0.05] <- "DOWN"
head(res_nLDL_vs_oxLDL3)
unique(res_nLDL_vs_oxLDL3$diffexpressed)
table(res_nLDL_vs_oxLDL3$diffexpressed)



#Create a new column "Top" that includes the names of the top differentially expressed genes

#ctrl vs nLDL
res_ctrl_vs_nLDL3$Top <- ifelse(res_ctrl_vs_nLDL3$row %in% head(arrange(res_ctrl_vs_nLDL3,pvalue),20)
                                $row & res_ctrl_vs_nLDL3$diffexpressed != "NO", res_ctrl_vs_nLDL3$row, NA)
head(res_ctrl_vs_nLDL3)

#ctrl vs oxLDL
res_ctrl_vs_oxLDL3$Top <- ifelse(res_ctrl_vs_oxLDL3$row %in% head(arrange(res_ctrl_vs_oxLDL3,pvalue),20)
                                 $row & res_ctrl_vs_oxLDL3$diffexpressed != "NO", res_ctrl_vs_oxLDL3$row, NA)
head(res_ctrl_vs_oxLDL3)

#nLDL vs oxLDL
res_nLDL_vs_oxLDL3$Top <- ifelse(res_nLDL_vs_oxLDL3$row %in% head(arrange(res_nLDL_vs_oxLDL3,pvalue),20)
                                 $row & res_nLDL_vs_oxLDL3$diffexpressed != "NO", res_nLDL_vs_oxLDL3$row, NA)
head(res_nLDL_vs_oxLDL3)

#plot

table(res_ctrl_vs_nLDL3$diffexpressed)
ctrl_vs_nLDL3 <- ggplot(res_ctrl_vs_nLDL3, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label=Top))+
  geom_vline(xintercept = 0, col = "grey", linetype = 'dashed', linewidth=0.8) +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = 'dashed', linewidth=0.8) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("aquamarine4", "grey", "sienna3"),
                     labels = c("DOWN", "NO", "UP")) +
  labs(x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"))+
  ggtitle('Ctrl vs nLDL-NMVs') +geom_text_repel(max.overlaps = Inf, size=2.5, fontface="bold")+
  theme(legend.title = element_blank(), plot.title.position = "right", plot.title = element_text(hjust = 0.5))+theme_classic()


table(res_ctrl_vs_oxLDL3$diffexpressed)
ctrl_vs_oxLDL3 <- ggplot(res_ctrl_vs_oxLDL3, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label=Top))+
  geom_vline(xintercept = 0, col = "grey", linetype = 'dashed', linewidth=0.8) +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = 'dashed', linewidth=0.8) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("aquamarine4", "grey", "sienna3"),
                     labels = c("DOWN", "NO", "UP")) +
  labs(x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"))+
  ggtitle('Ctrl vs oxLDL-NMVs') +geom_text_repel(max.overlaps = Inf, size=2.5, fontface="bold")+
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+theme_classic()

table(res_nLDL_vs_oxLDL3$diffexpressed)
nLDL_vs_oxLDL3 <- ggplot(res_nLDL_vs_oxLDL3, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label=Top))+
  geom_vline(xintercept = 0, col = "grey", linetype = 'dashed', linewidth=0.8) +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = 'dashed', linewidth=0.8) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("aquamarine4", "grey", "sienna3"),
                     labels = c("DOWN", "NO", "UP")) +
  labs(x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"))+
  ggtitle('nLDL-NMVs vs oxLDL-NMVs') +geom_text_repel(max.overlaps = Inf, size=2.5, fontface="bold")+
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+theme_classic()

ggarrange(ctrl_vs_nLDL3, ctrl_vs_oxLDL3,nLDL_vs_oxLDL3, nrow=1, common.legend = TRUE) %>% 
  annotate_figure(top = text_grob("Model: Condition + Participant",
                                  
                                                         face = "bold", size = 14))
# to Save it:
ggsave("Volcano plot condition plus participant.tiff", width = 15, height=6)


#========================================================================================


#Enrichment analysis

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")



library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(tidyverse)
library(ggpubr)



genes_ctrl_vs_nLDL3 <- sig_res_ctrl_vs_nLDL3[,1]
head(genes_ctrl_vs_nLDL3,20)

GO_BP_CTRL_vs_nLDL <- enrichGO (gene = genes_ctrl_vs_nLDL3, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO_CC_CTRL_vs_nLDL <- enrichGO (gene = genes_ctrl_vs_nLDL3, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "CC")
GO_MF_CTRL_vs_nLDL <- enrichGO (gene = genes_ctrl_vs_nLDL3, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "MF")

GO_BP_CTRL_vs_nLDL_df <- as.data.frame(GO_BP_CTRL_vs_nLDL)
GO_CC_CTRL_vs_nLDL_df <- as.data.frame(GO_CC_CTRL_vs_nLDL)
GO_MF_CTRL_vs_nLDL_df <- as.data.frame(GO_MF_CTRL_vs_nLDL)


BP_plot_CTRL_vs_nLDL <- head(GO_BP_CTRL_vs_nLDL_df,10) %>% ggplot(aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+
  coord_flip()+ggtitle("GO: Biological Process")
CC_plot_CTRL_vs_nLDL <- head(GO_CC_CTRL_vs_nLDL_df,10) %>% ggplot(aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+
  coord_flip()+ggtitle("GO: Cellular Component")
MF_plot_CTRL_vs_nLDL <- head(GO_MF_CTRL_vs_nLDL_df,10) %>% ggplot(aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+
  coord_flip()+ggtitle("GO: Molecular Function")

ggarrange(BP_plot, CC_plot,MF_plot, ncol=1, common.legend = TRUE) %>%
  annotate_figure(top = text_grob("Control vs nLDL MMVs" ,
                                  face = "bold", size = 14))
#my dot plots of GO


# Plot Biological Process (BP)
dotplot(GO_BP, title = "Biological Process", showCategory = 20)

# Plot Cellular Component (CC)
dotplot(GO_CC, title = "Cellular Component", showCategory = 20)

# Plot Molecular Function (MF)
dotplot(GO_MF, title = "Molecular Function", showCategory = 20)

# Combine the plots
install.packages("clusterProfiler")
install.packages("cowplot")

library(clusterProfiler)
library(cowplot)
combined_plot <- plot_grid(BP_plot, CC_plot, MF_plot, ncol = 1)

# Display the combined plot
print(combined_plot)




#KEGG

entrez_genes <- mapIds(org.Hs.eg.db, genes_ctrl_vs_nLDL3, 'ENTREZID', 'SYMBOL')
kegg_ctrl_vs_nLDL3 <- enrichKEGG(gene = entrez_genes, organism = 'hsa', keyType="kegg", pvalueCutoff = 0.05)
kegg_genes_ctrl_vs_nLDL_df <- as.data.frame(kegg_ctrl_vs_nLDL3)

#Plot KEGG
dotplot(kegg_ctrl_vs_nLDL3, showCategory=20)

#====================================================================================

# CTRL vs oxLDL Erich.

genes_ctrl_vs_oxLDL3 <- sig_res_ctrl_vs_oxLDL3[,1]
head(genes_ctrl_vs_oxLDL3,20)

GO_BP_CTRL_vs_oxLDL <- enrichGO (gene = genes_ctrl_vs_oxLDL3, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO_CC_CTRL_vs_oxLDL <- enrichGO (gene = genes_ctrl_vs_oxLDL3, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "CC")
GO_MF_CTRL_vs_oxLDL <- enrichGO (gene = genes_ctrl_vs_oxLDL3, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "MF")

GO_BP_CTRL_vs_oxLDL_df <- as.data.frame(GO_BP_CTRL_vs_oxLDL)
GO_CC_CTRL_vs_oxLDL_df <- as.data.frame(GO_CC_CTRL_vs_oxLDL)
GO_MF_CTRL_vs_oxLDL_df <- as.data.frame(GO_MF_CTRL_vs_oxLDL)


BP_plot_CTRL_vs_oxLDL <- head(GO_BP_CTRL_vs_oxLDL_df,10) %>% ggplot(aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+
  coord_flip()+ggtitle("GO: Biological Process")
CC_plot_CTRL_vs_oxLDL <- head(GO_CC_CTRL_vs_oxLDL_df,10) %>% ggplot(aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+
  coord_flip()+ggtitle("GO: Cellular Component")
MF_plot_CTRL_vs_oxLDL <- head(GO_MF_CTRL_vs_oxLDL_df,10) %>% ggplot(aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+
  coord_flip()+ggtitle("GO: Molecular Function")

ggarrange(BP_plot, CC_plot,MF_plot, ncol=1, common.legend = TRUE) %>%
  annotate_figure(top = text_grob("Control vs nLDL MMVs" ,
                                  face = "bold", size = 14))
#my dot plots of GO


# Plot Biological Process (BP)
dotplot(GO_BP, title = "Biological Process", showCategory = 20)

# Plot Cellular Component (CC)
dotplot(GO_CC, title = "Cellular Component", showCategory = 20)

# Plot Molecular Function (MF)
dotplot(GO_MF, title = "Molecular Function", showCategory = 20)

# Combine the plots
install.packages("clusterProfiler")
install.packages("cowplot")

library(clusterProfiler)
library(cowplot)
combined_plot <- plot_grid(BP_plot, CC_plot, MF_plot, ncol = 1)

# Display the combined plot
print(combined_plot)




#KEGG

entrez_genes <- mapIds(org.Hs.eg.db, genes_ctrl_vs_oxLDL3, 'ENTREZID', 'SYMBOL')
kegg_CTRL_vs_oxLDL3 <- enrichKEGG(gene = entrez_genes, organism = 'hsa', keyType="kegg", pvalueCutoff = 0.05)
kegg_kegg_CTRL_vs_oxLDL_df <- as.data.frame(kegg_CTRL_vs_oxLDL3)

#Plot KEGG
dotplot(kegg_CTRL_vs_oxLDL3, showCategory=20)


#================================================================================

# nLDL vs oxLDL Enrich.

genes_nLDL_vs_oxLDL3 <- sig_res_nLDL_vs_oxLDL3[,1]
head(genes_nLDL_vs_oxLDL3,20)

GO_BP_nLDL_vs_oxLDL <- enrichGO (gene = genes_nLDL_vs_oxLDL3, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO_CC_nLDL_vs_oxLDL <- enrichGO (gene = genes_nLDL_vs_oxLDL3, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "CC")
GO_MF_nLDL_vs_oxLDL <- enrichGO (gene = genes_nLDL_vs_oxLDL3, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "MF")

GO_BP_nLDL_vs_oxLDL_df <- as.data.frame(GO_BP_nLDL_vs_oxLDL)
GO_CC_nLDL_vs_oxLDL_df <- as.data.frame(GO_CC_nLDL_vs_oxLDL)
GO_MF_nLDL_vs_oxLDL_df <- as.data.frame(GO_MF_nLDL_vs_oxLDL)


BP_plot_nLDL_vs_oxLDL <- head(GO_BP_nLDL_vs_oxLDL_df,10) %>% ggplot(aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+
  coord_flip()+ggtitle("GO: Biological Process")
CC_plot_nLDL_vs_oxLDL <- head(GO_CC_nLDL_vs_oxLDL_df,10) %>% ggplot(aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+
  coord_flip()+ggtitle("GO: Cellular Component")
MF_plot_nLDL_vs_oxLDL <- head(GO_MF_nLDL_vs_oxLDL_df,10) %>% ggplot(aes(x=Description, y=Count, fill=p.adjust))+geom_bar(stat="identity")+
  coord_flip()+ggtitle("GO: Molecular Function")

ggarrange(BP_plot, CC_plot,MF_plot, ncol=1, common.legend = TRUE) %>%
  annotate_figure(top = text_grob("Control vs nLDL MMVs" ,
                                  face = "bold", size = 14))
#my dot plots of GO


# Plot Biological Process (BP)
dotplot(GO_BP, title = "Biological Process", showCategory = 20)

# Plot Cellular Component (CC)
dotplot(GO_CC, title = "Cellular Component", showCategory = 20)

# Plot Molecular Function (MF)
dotplot(GO_MF, title = "Molecular Function", showCategory = 20)

# Combine the plots
install.packages("clusterProfiler")
install.packages("cowplot")

library(clusterProfiler)
library(cowplot)
combined_plot <- plot_grid(BP_plot, CC_plot, MF_plot, ncol = 1)

# Display the combined plot
print(combined_plot)




#KEGG

entrez_genes <- mapIds(org.Hs.eg.db, genes_nLDL_vs_oxLDL3, 'ENTREZID', 'SYMBOL')
kegg_nLDL_vs_oxLDL3 <- enrichKEGG(gene = entrez_genes, organism = 'hsa', keyType="kegg", pvalueCutoff = 0.05)
kegg_nLDL_vs_oxLDL_df <- as.data.frame(kegg_nLDL_vs_oxLDL3)

#Plot KEGG
dotplot(kegg_nLDL_vs_oxLDL3, showCategory=20)

