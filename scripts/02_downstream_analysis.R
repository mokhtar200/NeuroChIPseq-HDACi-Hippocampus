# =============================================
# NeuroChIPseq HDACi - Downstream Analysis
# =============================================
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(DiffBind)
library(pheatmap)
library(rtracklayer)
library(GenomicRanges)
library(org.Mm.eg.db)
library(ggplot2)

# ------------------------------
# 1. Load Peaks
# ------------------------------
peak_file <- "../data/peaks/AcH2B_TSA_vs_Veh_peaks.broadPeak"
peaks <- import(peak_file)

# ------------------------------
# 2. Annotate Peaks
# ------------------------------
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peak_anno <- annotatePeak(peaks, TxDb=txdb, tssRegion=c(-3000,3000))
plotAnnoPie(peak_anno)

# ------------------------------
# 3. Prepare DiffBind analysis
# ------------------------------
samples <- data.frame(
  SampleID = c("Veh1","Veh2","TSA1"),
  Condition = c("Vehicle","Vehicle","TSA"),
  Replicate = c(1,2,1),
  Peaks = c("../data/peaks/Veh1_peaks.bed",
            "../data/peaks/Veh2_peaks.bed",
            "../data/peaks/TSA1_peaks.bed"),
  bamReads = c("../data/bam/SRR647745_sorted.bam",
               "../data/bam/SRR647746_sorted.bam",
               "../data/bam/SRR647747_sorted.bam"),
  stringsAsFactors = FALSE
)
dba_obj <- dba(sampleSheet = samples)
dba_obj <- dba.count(dba_obj)
dba_obj <- dba.contrast(dba_obj, categories=DBA_CONDITION, minMembers=1)
dba_obj <- dba.analyze(dba_obj)
db_results <- dba.report(dba_obj, method=DBA_DESEQ2, th=0.05)

# ------------------------------
# 4. Functional Enrichment
# ------------------------------
gene_ids <- as.data.frame(peak_anno)$geneId
ego <- enrichGO(gene=gene_ids, OrgDb=org.Mm.eg.db,
                keyType="ENTREZID", ont="BP",
                pAdjustMethod="BH", qvalueCutoff=0.05,
                readable=TRUE)
dotplot(ego, showCategory=20)

# ------------------------------
# 5. Heatmap of top peaks
# ------------------------------
top_peaks <- head(peak_anno, 50)
fold_changes <- as.data.frame(mcols(top_peaks))$Fold
pheatmap(matrix(fold_changes, ncol=1), cluster_rows=TRUE, cluster_cols=FALSE)

# ------------------------------
# 6. Save Results
# ------------------------------
write.csv(as.data.frame(db_results), file="../results/differential_binding/AcH2B_diff_binding.csv")
write.csv(as.data.frame(peak_anno), file="../results/peak_annotation/AcH2B_annotated_peaks.csv")
