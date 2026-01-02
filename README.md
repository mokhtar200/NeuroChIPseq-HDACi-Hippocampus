# NeuroChIPseq HDACi Hippocampus

This project analyzes the impact of **HDAC inhibition (TSA)** on **hippocampal histone acetylation** in mice using **ChIP-Seq**.  
We focus on the histone marks: **AcH2B, AcH3K9,14, AcH4K12**, comparing **TSA vs Vehicle** treatment.

## Dataset
- BioProject: [PRJNA186421](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA186421)
- SRA Accession IDs used: SRR647745-SRR647756 (TSA and Vehicle samples)
- Organism: Mus musculus
- Tissue: Hippocampus

## Workflow
1. Download raw SRR files.
2. Convert to FASTQ.
3. Align reads to **mm10** genome.
4. Call **broad peaks** using **MACS2**.
5. Annotate peaks and perform differential analysis in R.
6. Generate plots and functional enrichment results.

## Requirements

### Bash / Command-line:
- SRA Toolkit (`prefetch`, `fastq-dump`)
- bowtie2
- samtools
- MACS2

### R:
- ChIPseeker
- TxDb.Mmusculus.UCSC.mm10.knownGene
- clusterProfiler
- DiffBind
- pheatmap
- rtracklayer
- GenomicRanges

- `scripts/` : Bash and R scripts
- `results/` : annotated peaks, differential binding, heatmaps, and enrichment results
