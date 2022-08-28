#GSEA
library(tibble)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(patchwork)
library(GSVA)
library(pheatmap)
library(clusterProfiler)
library(stringr)
library(enrichplot)
library(org.Mm.eg.db)

setwd("D:\\Data\\scRNAseq\\scRNA PDL1 (LSEC)\\Tom+ GSEA")

data<-read.csv('./NASH_LSEC_DEG.csv')
#rank all DEGs
colnames(data)<-c("gene",'log2FC')
head(data)
ranks<- deframe(data)
ranks
##Setup MSigDB
msigdbr_species()
#create GMT file and DotPlot

m_df<- msigdbr(species = "Mus musculus", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
index <- lengths(fgsea_sets)>500
MmGMT <- fgsea_sets[!index]
MmGMT0 <- do.call(cbind, lapply(lapply(MmGMT, unlist), `length<-`, max(lengths(MmGMT))))
MmGMT0 <- as.matrix(t(MmGMT0))
MmGMT1 <- cbind(rownames(MmGMT0),"na",MmGMT0)
write.table(MmGMT1,"MmGMT1_H.gmt",sep = "\t",
            col.names = FALSE,row.names = FALSE,quote = FALSE)

m_df<- msigdbr(species = "Mus musculus", category = "C2", subcategory = 'CGP')
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
index <- lengths(fgsea_sets)>500
MmGMT <- fgsea_sets[!index]
MmGMT0 <- do.call(cbind, lapply(lapply(MmGMT, unlist), `length<-`, max(lengths(MmGMT))))
MmGMT0 <- as.matrix(t(MmGMT0))
MmGMT1 <- cbind(rownames(MmGMT0),"na",MmGMT0)
write.table(MmGMT1,"MmGMT1_C2_CGP.gmt",sep = "\t",
            col.names = FALSE,row.names = FALSE,quote = FALSE)


geneset <- read.gmt("MmGMT1_H.gmt")
egmt <- GSEA(ranks,pvalueCutoff = 0.15, TERM2GENE=geneset,verbose=F)
write.csv(egmt@result,file="GSEA_H.csv",row.names=TRUE)
saveRDS(egmt, file = "./Rdata/GSEA_H.rds")

geneset <- read.gmt("MmGMT1_C2_CGP.gmt") 
sene<-geneset[geneset$term=="FRIDMAN_SENESCENCE_UP",]
egmt <- GSEA(ranks,pvalueCutoff = 0.15, TERM2GENE=geneset,verbose=F)
write.csv(egmt@result,file="GSEA_C2_CGP.csv",row.names=TRUE)
saveRDS(egmt, file = "./Rdata/GSEA_C2_CGP.rds")

egmt<-readRDS( file = "./Rdata/GSEA_H.rds")
egmt<-readRDS( file = "./Rdata/GSEA_C2_CGP.rds")

View(egmt@result)

sene<-geneset[geneset$term=="FRIDMAN_SENESCENCE_UP",]
sene<-geneset[geneset$term=="FRIDMAN_SENESCENCE_DN",]
egmt <- GSEA(ranks,pvalueCutoff = 1, TERM2GENE=sene,verbose=F)

#SAVE dotplot of GSEA
tiff(file="GSEA_C7_ImmuneSigDB.tif",width = 1000, height = 500, res=100)
dotplot(egmt,split=".sign",orderBy = "p.adjust",font.size = 9,showCategory = 20)+facet_grid(~.sign)
dev.off()
x<-'HALLMARK_HYPOXIA'

#SAVE GSEA plot
tiff(file="HALLMARK_OXIDATIVE_PHOSPHORYLATION.tif",width = 2200, height = 1500, res=300)
gseaplot2(egmt, geneSetID = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
          title ="HALLMARK_OXIDATIVE_PHOSPHORYLATION", subplots = 1:3)
dev.off()
