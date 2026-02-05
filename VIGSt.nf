#!/usr/bin/env nextflow

// ===========================================
// VIGSt RNA-seq Pipeline
// ===========================================

nextflow.enable.dsl=2

// Define parameters
params.run_name = "Viroid_Viridiplantae_transcriptomic"
params.query = "(viroid[All Fields]) AND (Viridiplantae[Organism]) AND (transcriptomic[Source])"

//params.run_name = "VIGSt_Run_test"
//params.query = "SRR035182"
params.threads = 8


// 01 Get SRA IDs by querying SRA database
process _01_Query_SRA_IDs {
    cpus params.threads
    tag "$params.run_name"

    publishDir "$params.run_name/results/01_RNA_seq_raw_data", mode: "copy"

    input:
    val query

    output:
    path "viroid_RNA_Seq_runs_info.csv", emit: run_info
    path "viroid_RNA_Seq_samples_id.txt"

    script:
    """
    # Query SRA for run info
    esearch -db sra -query "$query" | efetch -format runinfo > viroid_RNA_Seq_runs_info.csv

    # Extract SRR IDs
    tail -n +2 viroid_RNA_Seq_runs_info.csv | cut -d',' -f1 > viroid_RNA_Seq_samples_id.txt    
    """
}

// 02 Download sample project metadata from SRA
process _02_Download_sample_project_metadata {
    cpus params.threads
    tag "$params.run_name"

    publishDir "$params.run_name/results/01_RNA_seq_raw_data", mode: "copy"

    input:
    val query

    output:
    path "viroid_RNA_Seq_metadata_project.csv"
    path "viroid_RNA_Seq_metadata_sample.csv"

    script:
    """
    # Download metadata for projects
    esearch -db sra -query "$query" \
        | elink -target bioproject \
        | efetch -format docsum \
        | xtract -pattern DocumentSummary -element Project_Acc,Project_Title,Project_Description \
        > viroid_RNA_Seq_metadata_project.csv

    # Download metadata for samples
    esearch -db sra -query "$query" \
        | elink -target biosample \
        | efetch -format docsum \
        | xtract -pattern DocumentSummary -element Accession,Title,Organism,Description \
        > viroid_RNA_Seq_metadata_sample.csv
    """
}

// 03 Download raw fastq files from SRA and separate into read size and each bioproject folder
process _03_Download_raw_fastq_files {
    cpus params.threads
    tag "$bioproject-$sample_ids"

    publishDir "$params.run_name/results/01_RNA_seq_raw_data/raw_samples", mode: "copy"

    input:
    tuple val(read_size), val(bioproject), val(sample_ids)

    output:
    path "${read_size}/${bioproject}/${sample_ids}*.fastq.gz", emit: raw_fastq_files

    script:
    """
    # Download raw fastq files using fasterq-dump
    fasterq-dump $sample_ids -e $task.cpus --outdir $read_size/$bioproject -p --split-files

    # Compress all FASTQ files in parallel
    pigz -p $task.cpus $read_size/$bioproject/$sample_ids*.fastq
    """
}

// 04 FastQC raw samples
process _04_FastQC_raw_samples {
    cpus params.threads
    tag "$bioproject-$fastq_short"

    publishDir "$params.run_name/results/02_quality_control/raw_data_short_reads", mode: "copy"

    input:
    tuple val(bioproject), path(fastq_short)

    output:
    path "$bioproject/*_fastqc.zip", emit: fastqc_zips
    path "$bioproject/*_fastqc.html"


    script:
    """
    mkdir -p $bioproject

    fastqc --threads $task.cpus --outdir $bioproject $fastq_short
    """
}

// 05 MultiQC raw samples
process _05_MultiQC_raw_samples {
    cpus params.threads
    tag "$bioproject"

    publishDir "$params.run_name/results/02_quality_control/raw_data_short_reads", mode: "copy"

    input:
    tuple val(bioproject), path(fastqc_zip_file)

    output:
    path "$bioproject/multiqc_report.html"

    script:
    """
    multiqc -o $bioproject $fastqc_zip_file
    """
}   

// 06 Nanoplot raw samples
process _06_NanoPlot_raw_samples {
    cpus params.threads
    tag "$bioproject-$sample_ids"

    publishDir "$params.run_name/results/02_quality_control/raw_data_long_reads", mode: "copy"

    input:
    tuple val(bioproject), val(sample_ids), path(fastq_long)

    output:
    path "${bioproject}/${sample_ids}_nanoplot_report"

    script:
    """
    NanoPlot --fastq $fastq_long -o ${bioproject}/${sample_ids}_nanoplot_report --threads $task.cpus
    """
}

// 07 fastp filtering for short reads
process _07_fastp_filtering {

    cpus params.threads
    tag "$bioproject-$sample_ids"

    publishDir "$params.run_name/results/03_filtered_reads/filtered_short_reads", mode: "copy"

    input:
    tuple val(bioproject), val(sample_ids), path(fastq_short_reads)

    output:
    path "${bioproject}/${sample_ids}*.fastq.gz", emit: filtered_short_fastq_files
    path "${bioproject}/${sample_ids}_fastp.html"

    script:
    """
    mkdir -p $bioproject

    if [ ${fastq_short_reads.size()} -eq 2 ]; then
        fastp \
            --thread $task.cpus \
            --in1 ${fastq_short_reads[0]} \
            --in2 ${fastq_short_reads[1]} \
            --out1 ${bioproject}/${sample_ids}_1.fastq.gz \
            --out2 ${bioproject}/${sample_ids}_2.fastq.gz \
            --trim_poly_g --trim_poly_x \
            --cut_front --cut_tail \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --length_required 10 \
            --html ${bioproject}/${sample_ids}_fastp.html 
    else
        fastp \
            --thread $task.cpus \
            --in1 ${fastq_short_reads[0]} \
            --out1 ${bioproject}/${sample_ids}.fastq.gz \
            --trim_poly_g --trim_poly_x \
            --cut_front --cut_tail \
            --qualified_quality_phred 20 \
            --length_required 10 \
            --html ${bioproject}/${sample_ids}_fastp.html 
    fi
    """
}

// 08 fastplong filtering for long reads
process _08_fastplong_filtering {

    cpus params.threads
    tag "$bioproject-$sample_ids"

    publishDir "$params.run_name/results/03_filtered_reads/filtered_long_reads", mode: "copy"

    input:
    tuple val(bioproject), val(sample_ids), path(fastq_long_reads)

    output:
    path "${bioproject}/${sample_ids}.fastq.gz", emit: filtered_long_fastq_files
    path "${bioproject}/${sample_ids}_fastp.html"

    script:
    """
    mkdir -p $bioproject
    
    fastplong \
        --thread $task.cpus \
        --in ${fastq_long_reads} \
        --out ${bioproject}/${sample_ids}.fastq.gz \
        --trim_poly_x \
        --cut_front --cut_tail \
        --cut_mean_quality 10 \
        --qualified_quality_phred 10 \
        --length_required 100 \
        --html ${bioproject}/${sample_ids}_fastp.html
    """
}


// Workflow block
workflow {
    _01_Query_SRA_IDs(params.query)

    _02_Download_sample_project_metadata(params.query)

    sample_meta_ch = _01_Query_SRA_IDs.out.run_info
        .splitCsv(header:true)
        .map {row -> 

        def sample_ids = row.Run.trim()
        def bioproject = row.BioProject.trim()
        def platform = row.Platform.trim()

        def read_size
        if (platform in ["ILLUMINA", "DNBSEQ", "BGISEQ", "ION_TORRENT"]) {
            read_size = "Short_reads"
        } else if (platform in ["OXFORD_NANOPORE", "PACBIO_SMRT"]) {
            read_size = "Long_reads"
        }

        tuple(read_size, bioproject, sample_ids)}

    _03_Download_raw_fastq_files(sample_meta_ch)   
    
    short_fastqc_ch = _03_Download_raw_fastq_files.out.raw_fastq_files
        .filter {fastq_short_file -> fastq_short_file.toString().contains("/Short_reads/")}
        .flatten()
        .map {fastq_short ->
            def bioproject = fastq_short.parent.name
            tuple(bioproject, fastq_short)}

    _04_FastQC_raw_samples(short_fastqc_ch)

    multiqc_ch = _04_FastQC_raw_samples.out.fastqc_zips
        .map {fastqc_zip_file -> 
            def bioproject = fastqc_zip_file.parent.name
            tuple(bioproject, fastqc_zip_file)}
            .groupTuple()

    _05_MultiQC_raw_samples(multiqc_ch)

    long_nanoplot_ch = _03_Download_raw_fastq_files.out.raw_fastq_files
        .filter {fastq_long_file -> fastq_long_file.toString().contains("/Long_reads/")}
        .flatten()
        .map {fastq_long ->
            def bioproject = fastq_long.parent.name
            def sample_ids = fastq_long.simpleName
            tuple(bioproject, sample_ids, fastq_long)}

    _06_NanoPlot_raw_samples(long_nanoplot_ch)

     filtered_short_fastq_ch = short_fastqc_ch
        .map {_bioproject, fastq_short ->
            def bioproject = fastq_short.parent.name
            def sample_ids = fastq_short.simpleName.replaceAll(/(_[12])$/, '')
            def fastq_short_reads = fastq_short
            tuple(bioproject, sample_ids, fastq_short_reads)}
        .groupTuple(by: [0,1])

    _07_fastp_filtering(filtered_short_fastq_ch)

    filtered_long_fastq_ch = long_nanoplot_ch
    
    _08_fastplong_filtering(filtered_long_fastq_ch)    
}

