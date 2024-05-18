library(ChIPseeker)
setwd("d:/")
peak_Tumor_vs_Spl <- readPeakFile("Tumor_vs_spl_F0_P0.05_forMerge.bed")
peak_no_diff <- readPeakFile("Merge_nodiff.bed")
peak_Spl_vs_Tumor <- readPeakFile("spl_vs_Tumor_F0_P0.05_forMerge.bed")

covplot(peak_Tumor_vs_Spl, chr = c("chr1", "chr2"))
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
peakAnno_Tumor_vs_Spl <- annotatePeak(
  peak_Tumor_vs_Spl,
  tssRegion = c(-2500, 500),
  TxDb =TxDb.Mmusculus.UCSC.mm10.knownGene,
  annoDb = "org.Mm.eg.db")

peakAnno_no_diffl <- annotatePeak(
  peak_no_diff,
  tssRegion = c(-2500, 500),
  TxDb =TxDb.Mmusculus.UCSC.mm10.knownGene,
  annoDb = "org.Mm.eg.db")

peakAnno_Spl_vs_Tumor <- annotatePeak(
  peak_Spl_vs_Tumor,
  tssRegion = c(-2500, 500),
  TxDb =TxDb.Mmusculus.UCSC.mm10.knownGene,
  annoDb = "org.Mm.eg.db")
plotAnnoPie(peakAnno_Tumor_vs_Spl)
plotAnnoPie(peakAnno_no_diffl)
plotAnnoPie(peakAnno_Spl_vs_Tumor)
