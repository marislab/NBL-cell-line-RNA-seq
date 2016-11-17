require(edgeR)
require(GenomicFeatures)
require(ggplot2)
require(reshape)

###########################################
##   from STAR output, generate counts   ##
###########################################

#set working directory
setwd("")

#list files with count data
files <- grep ("ReadsPerGene.out.tab",  dir("aligned_full"),value=TRUE)

# read each file containing count data into a list
cntlist<-list()
for(i in files){
  tab <-read.delim(paste("aligned_full/",i,sep=""),header=FALSE)
  ##used illumina tru-seq, so sense strand = column 4
  cntlist[[i]] <- tab[,4]
}
# colapse list of count vectors into matrix
cnt <- do.call(cbind,cntlist)
# change colnames and rownames
colnames(cnt) <- do.call(rbind,strsplit(colnames(cnt),"\\."))[,1]
rownames(cnt) <- tab[,1]

# remove non gene counts
cntgenes <- cnt[grep("N_",rownames(cnt),invert=TRUE),]


###########################################
##      Calculate Gene Lengths           ##
###########################################

# read GTF file
GTFfile <- "../shared_resources/gtf/refSeq_hg19_2016-03-03.gtf"
txdb <- makeTxDbFromGFF(GTFfile,format="gtf")
# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})


###########################################
##          Calculate FPKM               ##
###########################################

# calculate library sizes from original count data
libSize <- apply(cnt,2,sum)
# reorder gene lengths in a vector  of ncol(gene count matrix)
geneLength <- unlist(exonic.gene.sizes)[rownames(cntgenes)]
# calculate fpkm
fpkm_matrix <- rpkm(cntgenes,lib.size=libSize,gene.length=geneLength)

###########################################
##      Save count, fpkm files           ##
###########################################
#counts
write.table(cntgenes, paste(Sys.Date(), "-CellLineSTAR-counts_2pass_matrix.txt",
                            sep = ""), sep = "\t", quote = F)
#fpkm 
write.table(fpkm_matrix, paste(Sys.Date(), "-CellLineSTAR-fpkm_2pass_matrix.txt",
                                sep = ""), sep = "\t", quote = F)
#can save matrix as r data file
save(fpkm_matrix,file=paste(Sys.Date(), "-fpkm_matrix_all_celllines.rda", 
                             sep = ""))


###########################################
##DESEq2 Differential Expression Analysis##
###########################################

#read in MYCN status
mycn <- read.delim("2016-11-17-CellLine-MYCN-status.txt", sep = "\t", header = T)

##construct DESeq dataset for MYCN status:
dataset <- DESeqDataSetFromMatrix(countData = cntgenes,
                                  colData = mycn,
                                  design = ~ Status)

##remove rows with 0 reads
dataset <- dataset[ rowSums(counts(dataset)) > 0, ]

#amplified vs non-amp written 'non, amp' means +FC = greater in amp subset
dataset$Status <- factor(dataset$Status, levels=c("Nonamplified", "Amplified"))

#differential expression analysis
dataset <- DESeq(dataset)
res <- results(dataset)
resOrdered <- res[order(res$padj),]

#summary for distribution of DE-genes
summary(resOrdered)

#how many adj p-values were < 0.1?
sum(res$padj < 0.1, na.rm=TRUE)

#plot MA
plotMA(res, main="DESeq2", ylim=c(-2,2))

#reorder by  significance
res.1 <- subset(resOrdered, padj < 0.1)

#save table
write.table(as.data.frame(res.1),
            file="DiffExp_MYCN_adjp.10.txt", row.names = T, quote = F, sep = "\t")

#for list of genes
diffexpgenes <- sort(as.factor(rownames(res.1)))
