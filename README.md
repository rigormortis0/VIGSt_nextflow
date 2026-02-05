mamba env export --from-history > VIGSt_nextflow.yml
mamba env export > VIGSt_nextflow_full.yml


# SNAKEMAKE

mamba create -n VIGSt -c conda-forge -c nodefaults snakemake entrez-direct sra-tools bioconda::fastqc bioconda::multiqc bioconda::nanoplot

snakemake --cores 128

mamba activate VIGSt

snakemake --lint


# NEXTFLOW

mamba create -n VIGSt_nextflow -c conda-forge -c nodefaults bioconda::nextflow entrez-direct sra-tools bioconda::fastqc bioconda::multiqc bioconda::nanoplot bioconda::fastp bioconda::fastplong

mamba install bioconda::fastp

nextflow run VIGSt.nf -resume

nextflow lint VIGSt.nf

Todo list:
- Use Porechop for long reads
- Check short reads that became too short after removing primers, they are miRNA-Seq
- PRJNA646051 and PRJNA528793 has some control/infected samples in NCBI
