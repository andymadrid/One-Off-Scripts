# annotate dmrseq results
library(dmrseq)
library(ChIPseeker)
library(clusterProfiler)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
load("dmrseqResults.rdata")
sigRegions <- regions[which(regions$pval<0.05),]
dmrs2 <- sigRegions
peaks <- annotatePeak(dmrs2,tssRegion=c(-2000,2000),TxDb=txdb,annoDb="org.Mm.eg.db")
peaks
peaks <- as.data.frame(peaks)
all <- peaks[which(peaks$annotation != "Distal Intergenic"),]
e <- bitr(all$SYMBOL,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db")
genesAll <- e[,2]
egoAll <- enrichGO(gene=genesAll,ont="BP",readable=T,OrgDb="org.Rn.eg.db",pvalueCutoff=0.2,qvalueCutoff=1)
egoAll <- simplify(egoAll, cutoff=0.7, by="p.adjust", select_fun=min)
head(summary(egoAll))
