library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
txbd<-TxDb.Hsapiens.UCSC.hg38.knownGene
peak <- readPeakFile('D:/TEs_CoTEs_pipline/data/V2-TEs/non_olap/cores1000/85_Alu_Merged.bed')
covplot(peak)
peakAnno<-annotatePeak(peak,tssRegion=c(-3000,3000),TxDb=txbd,annoDb="org.Hs.eg.db")
plotAnnoPie(peakAnno)
peakAnno = as.data.frame(peakAnno)
