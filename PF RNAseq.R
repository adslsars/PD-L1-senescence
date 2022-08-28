which R
R --version
module avail R
module load R/4.0

R

quit()
as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo",intern=TRUE))

library(limma)
library(edgeR)
library(Glimma)
library(Rsubread)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ggplot2)
library(goseq)
library(fgsea)
library(NMF)
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(ggnewscale)
library(enrichplot)


setwd("/yshare1/ZETTAI_path_WA_slash_home_KARA/home/adslsars")
getwd()
fastq.files <- list.files(path = "./backup/RNAseq/PF RNA", pattern = ".fastq$", full.names = TRUE)
fastq.files

align(index="./RNAseq ref/mm10_tom/mm10_tom", 
      readfile1 = fastq.files[seq(1,11,2)], readfile2 = fastq.files[seq(2,12,2)],
      nthreads = 12)
bam.files <- list.files(path = "./backup/RNAseq/PF RNA", pattern = ".BAM$", full.names = TRUE)
bam.files
props <- propmapped(files=bam.files)
props

fc <- featureCounts(bam.files,isPairedEnd = TRUE, 
                    annot.ext="./RNAseq ref/mm10_tom/mm10.refGene.gtf/mm10.refGene.gtf",
                    isGTFAnnotationFile=TRUE,nthreads = 12)
save(fc,file="./backup/RNAseq/PF RNA/align_fc.Rdata")

fc.exon <- featureCounts(bam.files,isPairedEnd = TRUE, 
                         annot.ext="./RNAseq ref/mm10_tom/mm10.refGene.gtf/mm10.refGene.gtf",
                         allowMultiOverlap=TRUE,
                         isGTFAnnotationFile=TRUE,useMetaFeatures = FALSE,nthreads = 12)
save(fc.exon,file="./backup/RNAseq/PF RNA/align_fcexon.Rdata")

rm <- featureCounts(bam.files,isPairedEnd = TRUE, 
                         annot.ext="./RNAseq ref/mm10_tom/mm10.refGene.gtf/mm10_repeatMasker.gtf",
                         isGTFAnnotationFile=TRUE,nthreads = 12)
save(rm,file="./backup/RNAseq/PF RNA/align_rm.Rdata")



##Start
setwd("D:\\Data\\RNA seq data analysis\\hiseq\\PF RNAseq")
load("./align_rm.Rdata")
countdata <- rm$counts

load("./align_fcexon.Rdata")
countdata <- fc.exon$counts
sampleinfo<- data.frame( SampleName=c('Quie1','Quie2','Quie3','Sene1','Sene2','Sene3'),
                         CellType=c(rep('Quie',3),rep('Sene',3)),
                         Batch=gl(3, 1, 6,labels = c("1","2","3")))
sampleinfo
head(countdata)
colnames(countdata)
colnames(countdata) <- c('Quie1','Quie2','Quie3','Sene1','Sene2','Sene3')
table(colnames(countdata)==sampleinfo$SampleName)
countdata<-cbind(rownames(countdata),countdata)
p21<-countdata[countdata[,1]=="Cdkn1a",2:7]
barplot(p21[1])

load("./align_fc.Rdata")
countdata <- fc$counts
sampleinfo<- data.frame( SampleName=c('Quie1','Quie2','Quie3','Sene1','Sene2','Sene3'),
                         CellType=c(rep('Quie',3),rep('Sene',3)),
                         Batch=gl(3, 1, 6,labels = c("1","2","3")))
sampleinfo
head(countdata)
colnames(countdata)
colnames(countdata) <- c('Quie1','Quie2','Quie3','Sene1','Sene2','Sene3')
table(colnames(countdata)==sampleinfo$SampleName)# check colname equal to info
#Filter low express gene
myCPM <- cpm(countdata)
head(myCPM)
thresh <- myCPM > 0.1
head(thresh)
table(rowSums(thresh))
# at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
counts.keep <- countdata[keep,]
summary(keep)
dim(counts.keep)
par(mfrow=c(2,3))
plot(myCPM[,1],countdata[,1],xlab="CPM", ylab="Raw Count",ylim=c(0,200),xlim=c(0,3))
plot(myCPM[,2],countdata[,2],xlab="CPM", ylab="Raw Count",ylim=c(0,200),xlim=c(0,3))
plot(myCPM[,3],countdata[,3],xlab="CPM", ylab="Raw Count",ylim=c(0,200),xlim=c(0,3))
plot(myCPM[,4],countdata[,4],xlab="CPM", ylab="Raw Count",ylim=c(0,200),xlim=c(0,3))
plot(myCPM[,5],countdata[,5],xlab="CPM", ylab="Raw Count",ylim=c(0,200),xlim=c(0,3))
plot(myCPM[,6],countdata[,6],xlab="CPM", ylab="Raw Count",ylim=c(0,200),xlim=c(0,3))

#Convert counts to DGEList objec
dgeObj <- DGEList(counts.keep)
dgeObj
names(dgeObj)
dgeObj$samples
#Quality control
par(mfrow=c(1,1))
dgeObj$samples$lib.size
# The names argument tells the barplot to use the sample names on the x-axis
barplot(dgeObj$samples$lib.size, names=colnames(dgeObj), las=2)
title("Barplot of library sizes")
#SAVE
tiff(file="Barplot of library sizes.tif",width = 1000, height = 750,pointsize = 20)
barplot(dgeObj$samples$lib.size, names=colnames(dgeObj), las=2)
title("Barplot of library sizes")
dev.off()

# Get log2 counts per million
logcounts <- cpm(dgeObj,log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalized)")
#SAVE
tiff(file="Boxplots of logCPMs (unnormalized).tif",width = 1000, height = 750,pointsize = 20)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalized)")
dev.off()

#Multidimensional scaling plots
plotMDS(dgeObj)
par(mfrow=c(1,2))
levels(factor(sampleinfo$CellType))
col.cell <- c("purple","orange")[factor(sampleinfo$CellType)]
data.frame(sampleinfo$CellType,col.cell)

plotMDS(dgeObj,col=col.cell, pch=16,ylim=c(-4,4),xlim=c(-4,4))
legend("topright",cex=0.7,fill=c("purple","orange"),legend=levels(factor(sampleinfo$CellType)))
title("Cell type")
levels(factor(sampleinfo$Batch))
col.Batch <- c("blue","red","dark green")[factor(sampleinfo$Batch)]
col.Batch
plotMDS(dgeObj,col=col.Batch, pch=16,ylim=c(-4,4),xlim=c(-4,4))
legend("topright",cex=0.7,fill=c("blue","red","dark green"),legend=levels(factor(sampleinfo$Batch)))
title("Batch")

# Save the figure
tiff(file="MDSplot_Celltype.tif",width = 750, height = 750,pointsize = 20)
plotMDS(dgeObj,col=col.cell, pch=16,ylim=c(-4,4),xlim=c(-4,4))
legend("topright",fill=c("purple","orange"),legend=levels(factor(sampleinfo$CellType)))
title("Cell type")
dev.off()

#MDS with Glimma
labels <- paste(sampleinfo$SampleName, sampleinfo$CellType, sampleinfo$Batch)
group <- paste(sampleinfo$CellType,sampleinfo$Batch,sep=".")
group <- factor(group)
glMDSPlot(dgeObj, labels=labels, groups=group, folder="mds")

#Hierarchical clustering with heatmaps
var_genes <- apply(logcounts, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:200]
head(select_var)
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[sampleinfo$CellType]

# Plot the heatmap
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Top 500 most variable genes",
          ColSideColors=col.cell,scale="row")

# Save the heatmap
tiff(file="High_var_genes.heatmap.tif",width = 2000, height = 1500,pointsize = 20)
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes",ColSideColors=col.cell,scale="row")
dev.off()



#Normalisation for composition bias
# Apply normalisation to DGEList object
dgeObj <- calcNormFactors(dgeObj)
dgeObj$samples

par(mfrow=c(2,3))
plotMD(logcounts,column = 1)
abline(h=0,col="grey")
plotMD(logcounts,column = 3)
abline(h=0,col="grey")
plotMD(logcounts,column = 5)
abline(h=0,col="grey")
plotMD(logcounts,column = 2)
abline(h=0,col="grey")
plotMD(logcounts,column = 4)
abline(h=0,col="grey")
plotMD(logcounts,column = 6)
abline(h=0,col="grey")

par(mfrow=c(2,3))
plotMD(dgeObj,column = 1)
abline(h=0,col="grey")
plotMD(dgeObj,column = 3)
abline(h=0,col="grey")
plotMD(dgeObj,column = 5)
abline(h=0,col="grey")
plotMD(dgeObj,column = 2)
abline(h=0,col="grey")
plotMD(dgeObj,column = 4)
abline(h=0,col="grey")
plotMD(dgeObj,column = 6)
abline(h=0,col="grey")

# Save the figure
tiff(file="Distribution.tif",width = 1000, height = 750,pointsize = 20)
par(mfrow=c(2,3))
plotMD(dgeObj,column = 1)
abline(h=0,col="grey")
plotMD(dgeObj,column = 3)
abline(h=0,col="grey")
plotMD(dgeObj,column = 5)
abline(h=0,col="grey")
plotMD(dgeObj,column = 2)
abline(h=0,col="grey")
plotMD(dgeObj,column = 4)
abline(h=0,col="grey")
plotMD(dgeObj,column = 6)
abline(h=0,col="grey")
dev.off()

#Save data
save(group,dgeObj,sampleinfo,file="./fc_preprocessing.Rdata")

#Load data
load("./fc_preprocessing.Rdata")

# Create the two variables
group1 <- as.character(group)
type <- sapply(strsplit(group1, ".", fixed=T), function(x) x[1])
batch <- sapply(strsplit(group1, ".", fixed=T), function(x) x[2])

### Testing for differential expression
# 1. Fit the linear model
# Specify a design matrix with an intercept term
design.lm <- model.matrix(~ 0+type)
design.lm
### Voom transform the data
par(mfrow=c(1,1))
v <- voom(dgeObj,design.lm,plot = TRUE)

par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")

# Save the figure
tiff(file="Voom transform.tif",width = 2000, height = 1000,pointsize = 20)
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
abline(h=median(v$E),col="blue")
dev.off()


Search<-cbind(v$E, matrix(rownames(v$E)))


fit <- lmFit(v, design.lm)
names(fit)
head(coef(fit))

cont.matrix <- makeContrasts(SvsQ=typeSene - typeQuie,levels=design.lm)
cont.matrix
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
dim(fit.cont)
summa.fit <- decideTests(fit.cont, lfc=0.6)
summary(summa.fit)


glXYPlot(x=fit.cont$coefficients[,1], y=fit.cont$lods[,1],
         xlab="Coefficients", ylab="LOD score", main="SvsQ",
         counts=v$E, groups=group, status=summa.fit[,1],
         anno=fit.cont$genes, side.main="SYMBOL", folder="volcano-LM")

results.lm <- as.data.frame(fit.cont,n = Inf,row.names=rownames(fit.cont$coefficients))


ann <- select(org.Mm.eg.db,keytype="SYMBOL",keys=rownames(results.lm),columns=c("GENENAME","ENTREZID"))

table(ann$SYMBOL==rownames(results.lm))
results.lm <- cbind(ann,results.lm )

write.csv(results.lm,file="LM_DEG_Results_SvsQ.csv",row.names=TRUE)



#Export data
save(fit.cont,dgeObj,group,
     results.lm,summa.fit,
     file="./fc_DE.Rdata")
load("./fc_DE.Rdata")

data<-cbind(summa.fit,results.lm)

### GOseq analysis
up_genelist<-names(summa.fit@.Data[summa.fit@.Data[,1]==1,])
down_genelist<-names(summa.fit@.Data[summa.fit@.Data[,1]==-1,])
up_ego <- enrichGO(gene         = up_genelist,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05)
down_ego <- enrichGO(gene         = down_genelist,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = 'SYMBOL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff  = 0.05)

View(up_ego@result)
View(down_ego@result)

write.csv(up_ego@result,file="SvsQ_fc GOBP-up.csv",row.names=TRUE)
write.csv(down_ego@result,file="SvsQ_fc GOBP-down.csv",row.names=TRUE)

saveRDS(up_ego, file = "./SvsQ_fc_up_ego.rds")
saveRDS(down_ego, file = "./SvsQ_fc_down_ego.rds")

up_ego<-readRDS( file = "./SvsQ_fc_up_ego.rds")
down_ego<-readRDS( file = "./SvsQ_fc_down_ego.rds")

dotplot(up_ego,orderBy = "p.adjust",showCategory=20)
barplot(up_ego,showCategory=20)
heatplot(up_ego)
up_ego2 <- pairwise_termsim(up_ego)
emapplot(up_ego2,showCategory = 30)
emapplot_cluster(up_ego2,showCategory = 50,label_format = 20)


dotplot(down_ego,orderBy = "p.adjust",showCategory=20)
barplot(down_ego,showCategory=20)
heatplot(down_ego)
down_ego2 <- pairwise_termsim(down_ego)
emapplot(down_ego2,showCategory = 30)
emapplot_cluster(down_ego2,showCategory = 50,label_format = 20)


View(results.lm)

go.lm <- goana(fit.cont, coef="SvsQ",geneid = results.lm$ENTREZID,species = "Mm")
kegg.lm <- kegga(fit.cont, coef="SvsQ",geneid = results.lm$ENTREZID,species = "Mm")
write.csv(go.lm,file="LM_GO_Results_PvsN.csv",row.names=FALSE)
write.csv(kegg.lm,file="LM_KEGG_Results_PvsN.csv",row.names=FALSE)


# cneplot for interested GO terms (up-regulated DEG)
View(up_ego@result)
up_ego_int = up_ego[up_ego$ID %in% c(
  "GO:0071559", "GO:0006979","GO:0001819","GO:0030595","GO:0050863",
  "GO:0050727", "GO:0038061","GO:0032609","GO:0032606","GO:0090398",'GO:0002478'
), asis=T]

View(up_ego_int@result)


up_geneList_2 = results.lm$coefficients
names(up_geneList_2)=results.lm$SYMBOL
up_geneList_2 = sort(up_geneList_2,decreasing = T)

tiff(file="cnetplot-2.tif",width = 5000, height = 4000, res=300)
cnetplot(up_ego_int,showCategory = 11, categorySize="pvalue", foldChange=up_geneList_2,colorEdge = TRUE,
         cex_label_category=1.3, cex_label_gene=0.7 )+labs(title = "n-Sen GOBP")
dev.off()

cnetplot(up_ego_int, foldChange=up_geneList_2, circular = TRUE, 
         colorEdge = TRUE,node_label='gene',cex_label_gene = 0.7)
