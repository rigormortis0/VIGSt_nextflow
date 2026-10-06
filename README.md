# VIGSt

**Viroid Infection Gene-expression Study:** a Nextflow workflow for cross-study, cross-species meta-analysis of plant RNA-seq data from plant-viroid pathosystems.

## Overview

VIGSt integrates publicly available plant RNA-seq datasets to identify host genes, orthogroups, and biological processes associated with viroid infection.
The workflow performs study-specific differential-expression analysis and then compares responses across species at the orthogroup level.
Its outputs support selection of conserved and specific candidate host genes for downstream functional validation.

## Scientific Objectives

- Identify differentially expressed host genes within each plant-viroid study.
- Compare differentially regulated genes across species through orthogroups.
- Identify orthogroups with consistent upregulated or downregulated expression patterns.

## Workflow

The pipeline runs the following analysis stages:

1. Retrieve SRA run information and metadata, then apply study and sample exclusions.
2. Download reads and perform FastQC, MultiQC, NanoPlot, fastp, or fastplong quality control.
3. Map reads to viroid reference genomes with Bowtie2 and summarize viroid-mapping statistics.
4. Download or stage plant genomes and annotations, build STAR indices, and map reads to plant genomes.
5. Quantify gene-level reads with featureCounts and assemble count matrices by BioProject and species.
6. Create DESeq2 metadata from manually curated sample conditions and identify study-specific DEGs.
7. Extract proteins, infer orthogroups with OrthoFinder, and annotate proteins with eggNOG-mapper.
8. Summarize DEG-associated orthogroups across studies and identify conserved expression patterns.
9. Run TopGO enrichment and link significant GO terms back to orthogroups, genes, and eggNOG descriptions.


## Output Structure

The run directory contains the following main outputs:

- `01_RNA_seq_raw_data/`: SRA metadata, curated run information, and downloaded reads.
- `02_quality_control/`: FastQC, MultiQC, and NanoPlot reports.
- `03_filtered_reads/`: quality-filtered short reads.
- `04_viroids_mapping/`: viroid references, alignments, and mapping summaries.
- `05_plants_mapping/`: plant references, STAR indices, host alignments, and mapping statistics.
- `06_gene_count/`: featureCounts tables and assigned-read summaries.
- `07_differential_expression/`: count matrices, DESeq2 metadata, normalized counts, and DEG tables.
- `08_orthology/`: extracted proteins and OrthoFinder results.
- `09_eggnog/`: eggNOG annotations and DEG-associated orthogroup tables.
- `10_orthogroup_log2fc_summary/`: cross-study orthogroup expression summaries and pattern calls.
- `11_orthogroup_functional_enrichment/`: TopGO results and pattern-to-GO links.

## Key Result Files

- `07_differential_expression/deseq2_results/*_deseq2_results.csv`: study-specific differential-expression results.
- `09_eggnog/deseq2_eggnog_orthogroups/*_deseq2_eggnog_orthogroups.csv`: DEGs annotated with orthogroups, GO terms, KEGG terms, and descriptions.
- `10_orthogroup_log2fc_summary/orthogroup_cross_plant_gene_expression_patterns.csv`: orthogroups with consistent cross-study expression patterns.


