###########################################
## Correlate cell line and patient data  ##
###########################################
library(ggplot2)
setwd("")

#read in cell line and patient data, log 2 transform
celllines <- read.delim("2016-11-21-CellLineDiffExp_MYCN_adjp.10.txt", header = T,
                        sep = "\t")
celllines$gene_symbol <- rownames(celllines)
celllines$log2Exp <- log(celllines$baseMean, 2)

patients <- read.delim("~/Box Sync/Maris_Lab/Manuscripts/Harenza_NBLCellLines_ScientificData/Archive/2016-11-22-PtDiffExp_MYCN_adjp.10.txt", header = T,
                       sep = "\t")
patients$logExp <- log(patients$baseMean, 2)
patients$gene_symbol <- rownames(patients)

#merge files for correlation
both <- merge(celllines, patients, by="gene_symbol")

#pearson correlations
cor.test(both$log2Exp, both$logExp, method = "pearson")
cor.test(both$log2FoldChange.x, both$log2FoldChange.y, method = "pearson")

###########################################
##           Plot correlations           ##
###########################################

ggplot(both, aes(log2Exp, logExp)) + geom_point(alpha=.4) + 
  stat_smooth(method="lm", se=FALSE) +
  theme(panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(), 
    axis.text.y = element_text(colour="black", size = 18), 	  
    axis.text.x = element_text(colour="black", size = 18), 
    axis.title=element_text(colour="black", size = 24, face="bold")) +
  labs(y = expression("Patient log"[2]*" (BaseMean)"), x = expression("Cell Line log"[2]*" (BaseMean)"), title = "")

ggplot(both, aes(log2FoldChange.x, log2FoldChange.y)) + geom_point(alpha=.4) + 
  stat_smooth(method="lm", se=FALSE) +
  theme(panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(), 
    axis.text.y = element_text(colour="black", size = 18), 	  
    axis.text.x = element_text(colour="black", size = 18), 
    axis.title=element_text(colour="black", size = 24, face="bold")) +
  labs(y = expression("Patient log"[2]*" (Fold Change)"), x = expression("Cell Line log"[2]*" (Fold Change)"), title = "")

############################################
#Plot correlations of non-significant genes#
############################################

#read in/log 2 transform non-sign gene data from DESeq2 (filtered by adj p > 0.2)
cell.nonsig <- read.delim("CellLineDiffExp_MYCN_notDE2.txt", header = T,
                          sep = "\t")
cell.nonsig$gene_symbol <- rownames(cell.nonsig)
cell.nonsig$log2Exp <- log(cell.nonsig$baseMean, 2)

pt.nonsig <- read.delim("PtDiffExp_MYCN_notDE2.txt", header = T,
                        sep = "\t")
pt.nonsig$logExp <- log(pt.nonsig$baseMean, 2)
pt.nonsig$gene_symbol <- rownames(pt.nonsig)

#merge files
merge.nonsig <- merge(cell.nonsig, pt.nonsig, by="gene_symbol")

#pearson correlations
cor.test(merge.nonsig$log2Exp, merge.nonsig$logExp, method = "pearson")
cor.test(merge.nonsig$log2FoldChange.x, merge.nonsig$log2FoldChange.y, method = "pearson")

#plots
ggplot(merge.nonsig, aes(log2Exp, logExp)) + geom_point(alpha=.4) + 
  stat_smooth(method="lm", se=FALSE) +
  theme(panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(), 
    axis.text.y = element_text(colour="black", size = 18), 	  
    axis.text.x = element_text(colour="black", size = 18), 
    axis.title=element_text(colour="black", size = 24, face="bold")) +
  labs(y = expression("Patient log"[2]*" (BaseMean)"), x = expression("Cell Line log"[2]*" (BaseMean)"), title = "")

ggplot(merge.nonsig, aes(log2FoldChange.x, log2FoldChange.y)) + geom_point(alpha=.4) + 
  stat_smooth(method="lm", se=FALSE) +
  theme(panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(), 
    axis.text.y = element_text(colour="black", size = 18), 	  
    axis.text.x = element_text(colour="black", size = 18), 
    axis.title=element_text(colour="black", size = 24, face="bold")) +
  labs(y = expression("Patient log"[2]*" (Fold Change)"), x = expression("Cell Line log"[2]*" (Fold Change)"), title = "")
