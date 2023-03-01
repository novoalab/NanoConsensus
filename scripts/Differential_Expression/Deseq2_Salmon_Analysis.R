#################################################################################
####### Script to perform DE analysis by DeSeq2 using counts from Salmon ########
#################################################################################

#Import required libraries: 
library("stringr") 
library("DESeq2") 
library(tximportData)
library(tximport) 
library("ggplot2") 
library("EnhancedVolcano") 

#Read command-line arguments:
args = commandArgs(trailingOnly=TRUE)

#Read geneID table:
tx2gene <- read.table(args[2])  

#List samples and extract sample names:
ff <- list.files( path = ".", recursive=F, pattern = "*.counts$", full.names = TRUE )  
ff_names <- str_replace(ff, ".counts", "")  
ff_names <- str_replace(ff_names, "./", "") 
names(ff) <- ff_names 

#Create dataframe with condition and replicates data:  
pieces<-str_split_fixed(ff_names, "_",n=3) 
coldata<-data.frame("replicate" =pieces[, 3], "condition"=paste(pieces[, 1], pieces[, 2],sep="-")) 
row.names(coldata)<-colnames(txi.salmon$counts) 

#Extract condition names:
conditions <- unique(coldata$condition)

#Import salmon data:
txi.salmon <- tximport(ff, type = "salmon", tx2gene = tx2gene)

#Run DeSeq2: 
dds <- DESeqDataSetFromTximport(txi.salmon, colData=coldata, design=~condition) 
dds <- dds[ rowSums(counts(dds)) > 1, ] 
dds <- DESeq(dds) 

##Output the regressive log of the raw counts: 
rld <- rlog(dds, blind=FALSE) 
rld_data<-as.data.frame(assay(rld)) 
write.table(rld_data,file=paste("rlog_genes_", conditions[1],"_",conditions[2],".txt", sep=""), sep="\t", quote=FALSE)  

##Plotting PCA: 
pcaData <- plotPCA(rld, intgroup=c("replicate", "condition"), returnData=TRUE) 
percentVar <- round(100 * attr(pcaData, "percentVar")) 
pdf(paste("PCA_rlog_", conditions[1],"_",conditions[2],".pdf", sep=""), width=20, height=20) 
ggplot(pcaData, aes(PC1, PC2, color=condition, label=replicate)) + 
 geom_text(alpha = 0.8, size = 8, show.legend = FALSE) +  
 geom_point(size=2, alpha = 0.5) + 
 theme(legend.title=element_text(size=22), legend.text=element_text(size=20)) + 
 xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
 ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
 scale_x_reverse() 
dev.off() 


#Define accessory functions: 
ddsResult <- function(dds, contdef, cond1, cond2 ) { 
   filename<-paste(contdef, cond1, cond2, sep="_") 
   filename<-paste(filename, ".txt", sep="") 
   res<-results(dds, contrast=c(contdef,cond1,cond2)) 
   resOrdered <- res[order(res$padj),] 
   write.csv((resOrdered), file=filename, row.names=T, quote=F) 

   return (resOrdered) 

} 

#Generate text output and volcano plot:
res_C_N <- ddsResult(dds, "condition", conditions[1],conditions[2])

#Plotting:  
pdf(paste("volcano_", conditions[1],"_",conditions[2],".pdf", sep=""))
    print(EnhancedVolcano(res_C_N,
         title = paste(conditions[1]," versus ",conditions[2], sep=""), 
         lab = row.names(res_C_N), 
         x = 'log2FoldChange',  
         pCutoff = 0.05,  
         FCcutoff = 2,  
         drawConnectors = TRUE, 
         widthConnectors = 0.75, 
         maxoverlapsConnectors = Inf, 
         y = 'padj'))
dev.off() 
