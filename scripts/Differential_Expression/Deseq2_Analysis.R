########################################################
####### Script to perform DE analysis by DeSeq2 ########
########################################################

#Import required libraries: 
library("stringr") 
library("DESeq2") 
library("ggplot2") 
library("EnhancedVolcano") 

#Read command-line arguments:
args = commandArgs(trailingOnly=TRUE)

#Create count and conditions matrix: 
number<-2 
ff <- list.files( path = ".", recursive=F, pattern = args[1], full.names = TRUE ) 
counts.files <- lapply( ff, read.table) 
counts <- as.data.frame( sapply( counts.files, function(x) x[ , number ] ) ) 
fn <- basename(ff) 
fn <- gsub( ".counts", "", fn) 
colnames(counts) <- fn 
row.names(counts) <- counts.files[[1]]$V1 
counts2 <- counts[-grep("__",rownames(counts)),] 
pieces<-str_split_fixed(names(counts2), "_",n=3) 
coldata<-data.frame("replicate" =pieces[, 3], "condition"=paste(pieces[, 1], pieces[, 2],sep="-")) 
row.names(coldata)<-colnames(counts2) 

#Extract condition names:
conditions <- unique(coldata$condition)

#Run DeSeq2:
dds <- DESeqDataSetFromMatrix(countData = counts2, 
                             colData = coldata, 
                             design = ~ condition ) 
dds <- dds[ rowSums(counts(dds)) > 1, ] 
dds <- DESeq(dds) 

##Output the regressive log of the raw counts: 
rld <- rlog(dds, blind=FALSE) 
rld_data<-as.data.frame(assay(rld)) 

#If provided, load feature ID - feature name data:
if (length(args)>1){
    desc<-read.table(file=args[2], sep="\t", header=F)  
    desc2<-data.frame("ID"=desc[,1],"SEQ"=desc[ ,2]) 
    rld_genes<-merge(as.data.frame(rld_data),desc2, all=T, by.x=0, by.y="ID", sort=FALSE) 
    write.table(rld_genes, file=paste("rlog_genes_", conditions[1],"_",conditions[2],".txt", sep=""), sep="\t", quote=FALSE) 
 
} else {
    write.table(rld_data,file=paste("rlog_genes_", conditions[1],"_",conditions[2],".txt", sep=""), sep="\t", quote=FALSE)  
}

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

##Plotting volcano plots: 
#Depening on the availability of the feature ID - feature name data:
if (length(args)>1){
    ddsResult <- function(dds, desc, contdef, cond1, cond2 ) { 
        filename<-paste(contdef, cond1, cond2, sep="_") 
        filename<-paste(filename, ".txt", sep="") 
        res<-results(dds, contrast=c(contdef,cond1,cond2)) 
        resOrdered <- res[order(res$padj),]   
        resMerged<-merge(as.data.frame(resOrdered), desc2, by.x="row.names", by.y="ID", sort=FALSE) 
        resMerged$ID<-NULL 
        write.csv((resMerged),file=filename,row.names=F, quote=F) 
        return (resMerged) 

    } 
    
    res_C_N<-ddsResult(dds, desc2, "condition", conditions[1],conditions[2])
    volc_C_N<-res_C_N[, c("Row.names", "log2FoldChange", "pvalue", "padj","SEQ")] 

    pdf(paste("volcano_", conditions[1],"_",conditions[2],".pdf", sep=""))
    print(EnhancedVolcano(volc_C_N,  
        title = paste(conditions[1]," versus ",conditions[2], sep=""),  
        lab = volc_C_N$SEQ, 
        x = 'log2FoldChange',  
        pCutoff = 0.05,  
        FCcutoff = 1.5,  
        drawConnectors = TRUE, 
        widthConnectors = 0.75, 
        maxoverlapsConnectors = Inf, 
        y = 'padj'))
    dev.off()  
 
} else {
    ddsResult <- function(dds, contdef, cond1, cond2 ) { 
        filename<-paste(contdef, cond1, cond2, sep="_") 
        filename<-paste(filename, ".txt", sep="") 
        res<-results(dds, contrast=c(contdef,cond1,cond2)) 
        resOrdered <- res[order(res$padj),] 
        write.csv((resOrdered), file=filename, row.names=T, quote=F)       
        return (resOrdered) 

    }
    
    res_C_N <- ddsResult(dds, "condition", conditions[1],conditions[2])

    pdf(paste("volcano_", conditions[1],"_",conditions[2],".pdf", sep=""))
        print(EnhancedVolcano(res_C_N,
         title = paste(conditions[1]," versus ",conditions[2], sep=""), 
         lab = row.names(res_C_N), 
         x = 'log2FoldChange',  
         pCutoff = 0.05,  
         FCcutoff = 1.5,  
         drawConnectors = TRUE, 
         widthConnectors = 0.75, 
         maxoverlapsConnectors = Inf, 
         y = 'padj'))
    dev.off()  
}

