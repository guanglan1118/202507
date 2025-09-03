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




