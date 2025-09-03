# Step-by-step RNA-seq analysis pipeline (bulk)
## 0) Plan the experiment
- Replicates: ≥3 per condition is the sweet spot; avoid confounding (e.g., all cases on one day, all controls on another).

- Metadata table (colData) with columns like: sample_id, condition, batch, sex, RIN, library_type, pairing, strandedness.

- Reference: lock to one build/annotation (e.g., GRCh38 + GENCODE v4x); record exact versions.

## 1) Raw data QC
Goal: catch problems early.
- Run FastQC per FASTQ, aggregate with MultiQC.

- Key checks: per-base quality, adapter content, overrepresented sequences, duplication, insert size (PE), GC bias.

- Optional: FastQ Screen for contamination.

**CLI (bash)**
~~~
mkdir -p qc/fastqc
fastqc -t 8 -o qc/fastqc raw/*.fastq.gz
multiqc -o qc/multiqc qc/fastqc
~~~ 

## 2) Trimming 
- Use fastp or Trim Galore (adapters, low-quality tails). Skip if libraries are already clean.
~~~
fastp -i raw/S1_R1.fq.gz -I raw/S1_R2.fq.gz \
      -o clean/S1_R1.fq.gz -O clean/S1_R2.fq.gz \
      --detect_adapter_for_pe --thread 8 --html qc/fastp_S1.html
~~~

## 3) Quantification strategy (pick one)
A) Alignment-based (gene-level counts)

- Align with STAR (or HISAT2), then count with featureCounts.

- Pros: supports fusion/variant calling later; good for QC like junction saturation.

- Cons: slower/heavier.
~~~
# STAR index (once)
STAR --runThreadN 16 --runMode genomeGenerate \
     --genomeDir ref/star_GRCh38_gencodev4x \
     --genomeFastaFiles GRCh38.fa \
     --sjdbGTFfile gencode.v4x.annotation.gtf \
     --sjdbOverhang 100

# Align
STAR --runThreadN 16 --genomeDir ref/star_GRCh38_gencodev4x \
     --readFilesIn clean/S1_R1.fq.gz clean/S1_R2.fq.gz \
     --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts --outFileNamePrefix aln/S1_

# Count (exon→gene)
featureCounts -T 16 -p -s 0 \
  -a gencode.v4x.annotation.gtf -o counts/genes.txt aln/*Aligned.sortedByCoord.out.bam
~~~
**Post-alignment QC**: samtools stats/idxstats, RSeQC (e.g., infer_experiment.py for strandedness, geneBody_coverage.py), Picard (duplication).

## B) Alignment-free (transcript-level → gene)

- Quantify with Salmon (or kallisto), then collapse to genes via tximport in R.

- Pros: fast, accurate for DGE; handles strandedness autodetect (--libType A).

- Cons: less direct for fusions/variants.
~~~
# Packages
library(tximport); library(readr); library(DESeq2); library(ggplot2)
library(pheatmap); library(ComplexHeatmap); library(AnnotationDbi)
library(org.Hs.eg.db)  # or org.Mm.eg.db
# Metadata
coldata <- read.csv("metadata.csv")  # rows = samples; includes condition, batch, etc.

# ---- Branch A: featureCounts (gene x sample) ----
cts <- read.delim("counts/genes.txt", comment.char="#")
rownames(cts) <- cts$Geneid
cts <- cts[ , coldata$sample_id, drop=FALSE]  # order columns to match coldata

# ---- Branch B: Salmon (tximport) ----
files <- file.path("quant", coldata$sample_id, "quant.sf")
names(files) <- coldata$sample_id
tx2gene <- read.csv("tx2gene_gencodev4x.csv")  # transcript_id,gene_id
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

# Build DESeqDataSet (use txi$counts for Salmon)
dds <- DESeqDataSetFromMatrix(countData = if (exists("txi")) round(txi$counts) else cts,
                              colData = coldata,
                              design = ~ batch + condition)

# Prefilter low counts (speeds up)
keep <- rowSums(counts(dds) >= 10) >= max(2, floor(0.2 * ncol(dds)))
dds <- dds[keep,]

# Transform for QC
vsd <- vst(dds, blind=TRUE)

# PCA / outliers
plotPCA(vsd, intgroup=c("condition","batch"))
dists <- dist(t(assay(vsd)))
pheatmap(as.matrix(dists), clustering_distance_rows=dists, clustering_distance_cols=dists)

# Check library size, % mapped (add from MultiQC/STAR), strandedness notes
~~~
**Tips**:

- If UMI (3′ protocols), use a UMI-aware pipeline (e.g., STARsolo/alevin-fry) before tximport.

- If stranded, set featureCounts -s (0/1/2) correctly; Salmon -l A usually infers fine.

## 5) Model design & batch effects
- Put known covariates in the design: ~ batch + sex + RIN + condition.

- If unexplained structure remains, estimate surrogate variables (e.g., sva, RUVSeq) and include them.

- Avoid collinearity: do not include variables perfectly aligned with condition.

## 6) Differential expression (DE)
~~~
dds <- DESeq(dds)  # uses design above
res <- results(dds, contrast=c("condition","treated","control"))  # edit labels
res$symbol <- mapIds(org.Hs.eg.db, keys=row.names(res),
                     column="SYMBOL", keytype="ENSEMBL", multiVals="first")

# Shrink LFCs for ranking/plots
library(apeglm)
res_shr <- lfcShrink(dds, contrast=c("condition","treated","control"), type="apeglm")

# Filter/report
res_sig <- subset(as.data.frame(res_shr), padj < 0.05 & abs(log2FoldChange) >= 1)
write.csv(res_sig, "DE_treated_vs_control.csv", row.names=FALSE)
~~~
**Diagnostics**

- MA plot (plotMA(res_shr)), p-value histogram, mean-variance trend.

- Independent filtering is automatic in DESeq2; report the alpha used.

## 7) Visualization & summaries

- Volcano: log2FC vs −log10(padj).

- Heatmaps: top DE genes, z-scored VST expression with sample annotations (condition/batch).

- Pathways:

  - GSEA with fgsea using a ranked vector (e.g., stat = sign(log2FC)*-log10(pvalue) or lfcShrink).

  - ORA with clusterProfiler (GO/KEGG/Reactome/MSigDB).

## 8) Optional analyses (branch as needed)

- Isoform/DTU: DEXSeq, DRIMSeq, satuRn.

- Fusions: STAR-Fusion or Arriba (requires alignment BAMs).

- Variant calling: GATK RNA-seq best practices (split-NCigar, etc.).

- Deconvolution (bulk): MuSiC, CIBERSORTx, xCell.

- Gene set variation per sample: GSVA.

## 9) Reproducibility & packaging

- Workflow engine: Snakemake/Nextflow; version-pin with conda/mamba or Docker.

- In R, capture package versions with renv. Knit a Quarto/R Markdown report (methods, versions, figures).

- Save all: MultiQC, STAR logs, counting logs, sessionInfo(), and exact reference files.

- Consider GEO/SRA deposition; include raw FASTQs, processed counts, metadata, and a README.

## 10) Practical thresholds & gotchas

- Mapping rate: typically >70% (polyA, human). Much lower → check contamination/annotation mismatch.

- Strandedness: mis-specification can invert DE; verify with RSeQC infer_experiment.py.

- Low RIN/degradation: inspect gene body 5’–3’ bias; consider including RIN as a covariate.

- Low counts filtering: require expression in ≥20% of samples or ≥10 reads in ≥2 samples (tune to depth).

- Multiple contrasts: use a single fitted model; extract contrasts with results() per comparison, control FDR per contrast.

## 11) Minimal “starter” file list to archive

- metadata.csv, tx2gene_*.csv (if tximport), exact FASTA/GTF or Salmon index, counts matrix, vsd.rds, dds.rds, DE_*.csv, GSEA_results.csv, all QC HTML (FastQC, MultiQC), and an analysis.qmd with sessionInfo().














