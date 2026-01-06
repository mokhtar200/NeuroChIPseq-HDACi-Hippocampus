# ---------------------------
# 1️⃣ Load libraries
# ---------------------------
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(clusterProfiler)
library(GenomicRanges)
library(ggplot2)
library(VennDiagram)

# ---------------------------
# 2️⃣ Set working directories
# ---------------------------
# Change this path to where your peaks folders are
base_dir <- "~/Desktop/histone"
marks <- c("AcH2B","AcH3K9,14","AcH4K12")

# Initialize TxDb
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# ---------------------------
# 3️⃣ Function to analyze a histone mark
# ---------------------------
analyze_mark <- function(mark){
  
  cat("\n\n==== Processing:", mark, "====\n\n")
  
  peak_file <- file.path(base_dir, mark, "peaks", paste0(mark, "_TSA_vs_Veh_peaks.broadPeak"))
  
  # 1️⃣ Read peak file
  peaks <- readPeakFile(peak_file)
  
  # 2️⃣ Annotate peaks
  peak_anno <- annotatePeak(peaks,
                            TxDb = txdb,
                            tssRegion = c(-3000, 3000),
                            annoDb = "org.Mm.eg.db")
  
  # 3️⃣ Save annotation table
  out_table <- file.path(base_dir, mark, "R_analysis", paste0(mark, "_peak_annotation.csv"))
  dir.create(dirname(out_table), showWarnings = FALSE)
  write.csv(as.data.frame(peak_anno), out_table, row.names = FALSE)
  
  # 4️⃣ Plot annotation pie
  pie_plot_file <- file.path(base_dir, mark, "R_analysis", paste0(mark, "_annotation_pie.png"))
  png(pie_plot_file, width = 800, height = 600)
  plotAnnoPie(peak_anno)
  dev.off()
  
  # 5️⃣ Distance to TSS plot
  tss_plot_file <- file.path(base_dir, mark, "R_analysis", paste0(mark, "_TSS_distance.png"))
  png(tss_plot_file, width = 800, height = 600)
  plotDistToTSS(peak_anno, title=paste(mark, "Peaks around TSS"))
  dev.off()
  
  # 6️⃣ Extract gene IDs for GO enrichment
  gene_ids <- as.data.frame(peak_anno)$geneId
  
  # 7️⃣ GO enrichment (Biological Process)
  ego <- enrichGO(gene = gene_ids,
                  OrgDb = org.Mm.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05)
  
  # Save GO results
  go_file <- file.path(base_dir, mark, "R_analysis", paste0(mark, "_GO_BP.csv"))
  write.csv(as.data.frame(ego), go_file, row.names = FALSE)
  
  # Dotplot of top GO terms
  dotplot_file <- file.path(base_dir, mark, "R_analysis", paste0(mark, "_GO_BP_dotplot.png"))
  png(dotplot_file, width = 1000, height = 800)
  print(dotplot(ego, showCategory = 20))
  dev.off()
  
  return(list(peakAnno=peak_anno, GO=ego))
}

# ---------------------------
# 4️⃣ Run analysis for all marks
# ---------------------------
results <- lapply(marks, analyze_mark)

names(results) <- marks

# ---------------------------
# 5️⃣  Venn diagram of peak-associated genes between marks
# ---------------------------
# Extract gene lists
gene_lists <- lapply(results, function(x) as.data.frame(x$peakAnno)$geneId)

# Plot Venn diagram
venn_file <- file.path(base_dir, "venn_peak_genes.png")
png(venn_file, width = 800, height = 800)
venn.plot <- venn.diagram(
  x = gene_lists,
  category.names = marks,
  filename = NULL,
  fill = c("red","green","blue"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5
)
grid.draw(venn.plot)
dev.off()

cat("\n\n==== All analysis complete! Files saved under R_analysis folders ==== \n")
