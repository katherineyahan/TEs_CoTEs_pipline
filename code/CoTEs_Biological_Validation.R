library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(gprofiler2)
txdb<-TxDb.Hsapiens.UCSC.hg38.knownGene
peak <- readPeakFile('D:/TEs_CoTEs_pipline/data/V2-TEs/non_olap/cores1000/542_L1MA1_Merged.bed')
covplot(peak)
peakAnno<-annotatePeak(peak,tssRegion=c(-3000,3000),TxDb=txdb,annoDb="org.Hs.eg.db")
plotAnnoPie(peakAnno)
peakAnno = as.data.frame(peakAnno)

filtered_df<-subset(peakAnno,peakAnno$seqnames=="chrX")
genes_list<-filtered_df$SYMBOL
gostres<-gost(query=c(genes_list),organism="hsapiens",significant=TRUE,user_threshold = 0.05)
