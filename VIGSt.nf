#!/usr/bin/env nextflow

// ===========================================
// VIGStargets RNA-seq Pipeline
// ===========================================

nextflow.enable.dsl=2

// Define parameters
params.run_name = "Run_complete"
params.query = "(viroid[All Fields]) AND (Viridiplantae[Organism]) AND (transcriptomic[Source])"
params.threads = 8

// 01 Get SRA IDs by querying SRA database
process _01_Query_SRA_IDs {

    cpus params.threads
    tag "SRA IDs"

    publishDir "$params.run_name/01_RNA_seq_raw_data", mode: "copy"

    input:
    val query

    output:
    path "viroid_RNA_Seq_runs_info.csv", emit: run_info

    script:
    """
    # Query SRA for run info
    esearch -db sra -query "$query" | efetch -format runinfo > viroid_RNA_Seq_runs_info.csv
    """
}

// 02 Exclude samples from run 
process _02_Exclude_samples {

    cpus params.threads
    tag "Exclude samples"

    publishDir "$params.run_name/01_RNA_seq_raw_data", mode: "copy"

    input:
    path run_info
    path excluded_samples

    output:
    path "viroid_RNA_Seq_runs_info_clean.csv", emit: run_info_clean
    path "viroid_RNA_Seq_runs_info_simplified_clean.csv", emit: run_info_simplified_clean

    script:
    """
    awk -F',' '
        NR==FNR {
            if (FNR > 1) excluded_samples[\$1]
            next
        }
        FNR==1 {
            print
            next
        }
        !(\$1 in excluded_samples)
    ' ${excluded_samples} ${run_info} > viroid_RNA_Seq_runs_info_clean.csv
   
    # Keep only columns of interest 
    cut -d',' -f1,4,7,12,13,19,22,29,30 viroid_RNA_Seq_runs_info_clean.csv > viroid_RNA_Seq_runs_info_simplified_clean.csv   
    """
}

// 03 Exclude complete studies from run
process _03_Exclude_studies {

    cpus params.threads
    tag "Exclude studies"

    publishDir "$params.run_name/01_RNA_seq_raw_data", mode: "copy"

    input:
    path run_info
    path excluded_studies

    output:
    path "viroid_RNA_Seq_runs_info_studies_clean.csv", emit: run_info_studies_clean

    script:
    """
    awk -F',' '
        function trim(x) {
            gsub(/^[[:space:]]+|[[:space:]]+\$/, "", x)
            return x
        }
        NR==FNR {
            if (FNR > 1) excluded_studies[trim(\$1)]
            next
        }
        FNR==1 {
            print
            next
        }
        !(trim(\$22) in excluded_studies)
    ' ${excluded_studies} ${run_info} > viroid_RNA_Seq_runs_info_studies_clean.csv
    """
}

// 04 Download sample project metadata from SRA
process _04_Download_sample_project_metadata {

    cpus params.threads
    tag "Download metadata"

    publishDir "$params.run_name/01_RNA_seq_raw_data", mode: "copy"

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

// 05 Download raw fastq files from SRA and separate into read size and each bioproject folder
process _05_Download_raw_fastq_files {

    cpus params.threads
    tag "$bioproject-$sample_ids"

    publishDir "$params.run_name/01_RNA_seq_raw_data/raw_samples", mode: "copy"

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

// 06 FastQC raw samples
process _06_FastQC_raw_samples {

    cpus params.threads
    tag "$bioproject-$fastq_short"

    publishDir "$params.run_name/02_quality_control/raw_data_short_reads", mode: "copy"

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

// 07 MultiQC raw samples
process _07_MultiQC_raw_samples {

    cpus params.threads
    tag "$bioproject"

    publishDir "$params.run_name/02_quality_control/raw_data_short_reads", mode: "copy"

    input:
    tuple val(bioproject), path(fastqc_zip_file)

    output:
    path "$bioproject/multiqc_report.html"

    script:
    """
    multiqc -o $bioproject $fastqc_zip_file
    """
}   

// 08 Nanoplot raw samples
process _08_NanoPlot_raw_samples {

    cpus params.threads
    tag "$bioproject-$sample_ids"

    publishDir "$params.run_name/02_quality_control/raw_data_long_reads", mode: "copy"

    input:
    tuple val(bioproject), val(sample_ids), path(fastq_long)

    output:
    path "${bioproject}/${sample_ids}_nanoplot_report"

    script:
    """
    NanoPlot --fastq $fastq_long -o ${bioproject}/${sample_ids}_nanoplot_report --threads $task.cpus
    """
}

// 09 Fastp filtering for short reads
process _09_Fastp_filtering {

    cpus params.threads
    tag "$bioproject-$sample_ids"

    publishDir "$params.run_name/03_filtered_reads/filtered_short_reads", mode: "copy"

    input:
    tuple val(bioproject), val(sample_ids), path(fastq_short_reads)

    output:
    tuple val(bioproject), val(sample_ids), path("${bioproject}/${sample_ids}*.fastq.gz"), emit: filtered_short_fastq_files
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

// 10 Fastplong filtering for long reads
process _10_Fastplong_filtering {

    cpus params.threads
    tag "$bioproject-$sample_ids"

    publishDir "$params.run_name/03_filtered_reads/filtered_long_reads", mode: "copy"

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

// 11 Download viroid refseq genomes
process _11_Download_viroid_genomes {

    cpus params.threads
    tag "Viroid_genomes"

    publishDir "$params.run_name/04_viroids_mapping", mode: "copy"

    output:
    path "viroid_refseq.fasta", emit: viroid_fasta

    script:
    """
    esearch -db nucleotide -query "Viroids[Organism] AND srcdb_refseq[PROP] AND "complete genome"[Title]" \
    | efetch -format fasta \
    > viroid_refseq_ori.fasta

    sed '/^>/ s/ /_/g' viroid_refseq_ori.fasta > viroid_refseq.fasta
    """
}

// 12 Build bowtie2 index for viroid genomes
process _12_Build_viroid_index {

    cpus params.threads
    tag "Viroid_index"

    publishDir "$params.run_name/04_viroids_mapping/index", mode: "copy"

    input:
    path viroid_fasta

    output:
    path "viroid_index.*", emit: viroid_index

    script:
    """
    unset PERL5LIB
    unset PERLLIB

    bowtie2-build ${viroid_fasta} viroid_index
    """
}

// 13 Map samples into viroid genomes and detect their infection using bowtie2
process _13_Map_to_viroid_genomes {

    cpus params.threads
    tag "$bioproject-$sample_ids"

    publishDir "$params.run_name/04_viroids_mapping/mapped", mode: "copy"

    input:
    tuple val(bioproject), val(sample_ids), path(filtered_fastq_files)
    path viroid_index

    output:
    path "${sample_ids}_viroid.bam", emit: aligned_bam

    script:
    """
    unset PERL5LIB
    unset PERLLIB

    if [ ${filtered_fastq_files.size()} -eq 2 ]; then
        bowtie2 -x viroid_index \
            -1 ${filtered_fastq_files[0]} \
            -2 ${filtered_fastq_files[1]} \
            --very-sensitive \
            -p $task.cpus \
        | samtools sort -@ $task.cpus -o ${sample_ids}_viroid.bam
    else
        bowtie2 -x viroid_index \
            -U ${filtered_fastq_files[0]} \
            --very-sensitive \
            -p $task.cpus \
        | samtools sort -@ $task.cpus -o ${sample_ids}_viroid.bam
    fi

    samtools index ${sample_ids}_viroid.bam
    """
}

// 14 Summarize statistics of sample reads mapping into viroid genomes  
process _14_Map_to_viroid_genomes_stats {

    cpus params.threads
    tag "Map_statistics_summary"

    publishDir "$params.run_name/04_viroids_mapping/mapped", mode: "copy"

    input:
    tuple val(bioproject), val(sample_ids), path(filtered_fastq_files)
    path aligned_bam

    output:
    path "${sample_ids}_viroid_read_summary.csv", emit: virod_stats
    
    script:
    """
    # Total reads per sample
    if [ ${filtered_fastq_files.size()} -eq 2 ]; then
        r1=\$(zcat ${filtered_fastq_files[0]} | wc -l)
        r2=\$(zcat ${filtered_fastq_files[1]} | wc -l)
        total_reads=\$(( (r1 + r2) / 4 ))
    else
        r1=\$(zcat ${filtered_fastq_files[0]} | wc -l)
        total_reads=\$(( r1 / 4 ))
    fi

    # Total reads mapped to viroids per sample
     viroid_reads=\$(samtools view -c -F 4 ${aligned_bam})

    # Percentage of mapped reads
    percent=\$(awk -v v=\$viroid_reads -v t=\$total_reads 'BEGIN {printf "%.5f", (t>0 ? v/t*100 : 0)}')

    # Summary tables per read
    echo "Run,total_reads,viroid_mapped_reads,percent" > ${sample_ids}_viroid_read_summary.csv
    echo "${sample_ids},\$total_reads,\$viroid_reads,\$percent" >> ${sample_ids}_viroid_read_summary.csv
    """
}

// 15 Summarize statistics of sample reads mapping into viroid genomes  
process _15_Merge_viroid_csvs {

    cpus params.threads
    tag "Map_statistics_summary"
    
    publishDir "$params.run_name/04_viroids_mapping", mode: "copy"

    input:
    path virod_tables

    output:
    path "viroid_summary_merged.csv", emit: viroid_map_stats

    script:
    """
    # Take header from first file
    head -n 1 ${virod_tables[0]} > viroid_summary_merged.csv

    # Append all rows (skip headers)
    for f in ${virod_tables}; do
        tail -n +2 \$f >> viroid_summary_merged.csv
    done
    """
}

// 16 Download reference plant genomes and gtf annotations from NCBI 
process _16_Download_plant_genomes {

    cpus params.threads
    tag "$plant_species"

    publishDir "$params.run_name/05_plants_mapping/genomes", mode: "copy"

    input:
    tuple val(plant_species), val(genome_url), val(gtf_url)

    output:
    tuple val(plant_species), path("${plant_species}/genome.fna"), path("${plant_species}/annotation.gtf"), emit: genome_gtf

    script:
    """
    mkdir -p ${plant_species}

    wget "$genome_url"
    wget "$gtf_url"

    gunzip -c *genomic.fna.gz > ${plant_species}/genome.fna
    gunzip -c *genomic.gtf.gz > ${plant_species}/annotation.gtf
    """
}

// 17 Downloaded plant genomes from diferent sources 
process _17_Downloaded_plant_genomes{

    cpus params.threads
    tag "$plant_species"

    publishDir "$params.run_name/05_plants_mapping/genomes", mode: "copy"

    input:
    tuple val(plant_species), path(genome_fna), path(annotation_gff3)

    output:
    tuple val(plant_species), path("${plant_species}/genome.fna"), path("${plant_species}/annotation.gtf"), emit: genome_gtf

    script:
    """
    mkdir -p ${plant_species}
    
    cp ${genome_fna} ${plant_species}/genome.fna

    gffread ${annotation_gff3} -T -o ${plant_species}/annotation.gtf
    """
}

// 18 Build start index for plant genomes
process _18_Build_plant_index {

    cpus params.threads
    tag "$plant_species"

    publishDir "$params.run_name/05_plants_mapping/index", mode: "copy"

    input:
    tuple val(plant_species), path(genome_fasta_file), path(annotation_gtf_file)

    output:
    tuple val(plant_species), path("${plant_species}_index"), emit: plant_index
    
    script:
    """
    mkdir -p ${plant_species}

    STAR \
        --runThreadN $task.cpus \
        --runMode genomeGenerate \
        --genomeDir ${plant_species}_index \
        --genomeFastaFiles $genome_fasta_file \
        --sjdbGTFfile $annotation_gtf_file \
        --limitGenomeGenerateRAM 90000000000 \
        --sjdbOverhang 149
    """
}

// 19 Map samples into plant genomes using STAR
process _19_Map_to_plant_genomes {

    cpus params.threads
    tag "$bioproject-$sample_ids-$plant_species"

    publishDir "$params.run_name/05_plants_mapping/mapped", mode: "copy"

    input:
    tuple val(plant_species), val(bioproject), val(sample_ids), path(filtered_fastq_files), val(single_pair_end), path(plant_index)

    output:
    tuple val(plant_species), val(bioproject), val(sample_ids), val(single_pair_end), path("${sample_ids}_Aligned.sortedByCoord.out.bam"), emit: aligned_bam
    path "${sample_ids}_Log.final.out", emit: log_map_stats

    script:
    """
    if [ ${filtered_fastq_files.size()} -eq 2 ]; then
        STAR \
            --runThreadN $task.cpus \
            --genomeDir $plant_index \
            --readFilesIn ${filtered_fastq_files[0]} ${filtered_fastq_files[1]} \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${sample_ids}_ 
    else
        STAR \
            --runThreadN $task.cpus \
            --genomeDir $plant_index \
            --readFilesIn ${filtered_fastq_files[0]} \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${sample_ids}_ 
    fi
    """ 
}

// 20 Summarize statistics of sample reads mapping into plant genomes
process _20_Map_to_plant_genomes_stats {

    cpus params.threads
    tag "Map_statistics_summary"

    publishDir "$params.run_name/05_plants_mapping", mode: "copy"

    input:
    path log_map_stats

    output:
    path "Map_stats_summary.csv", emit: plant_map_stats

    script:
    """
    echo "Run,Number_filtered_reads,Uniquely_mapped_perc,Multi_mapped_perc,Unmapped_short_perc,Unmapped_other_perc" > Map_stats_summary.csv

    # Loop through all STAR log files
    for log in ${log_map_stats}; do

        sample=\$(basename \$log | sed 's/_Log.final.out//')

        awk -F '|' -v sample_ids=\$sample '
            BEGIN { OFS="," }

            /Number of input reads/ {total=\$2+0}
            /Uniquely mapped reads %/ {unique=\$2+0}
            /% of reads mapped to multiple loci/ {multi=\$2+0}
            /% of reads unmapped: too short/ {short=\$2+0}
            /% of reads unmapped: other/ {other=\$2+0}
            
            END {
                print sample_ids, total, unique, multi, short, other
            }
        ' \$log >> Map_stats_summary.csv

    done
    """ 
}

// 21 Count mapped reads using featureCounts
process _21_Read_count_mapped {

    cpus params.threads
    tag "$bioproject-$sample_ids-$plant_species"

    publishDir "$params.run_name/06_gene_count/read_counts", mode: "copy"

    input:
    tuple val(plant_species), val(bioproject), val(sample_ids), val(single_pair_end), path(aligned_bam), path(annotation_gtf_file)

    output:
    tuple val(plant_species), val(bioproject), val(sample_ids), path("${sample_ids}_featureCounts.txt"), emit: read_counts
    path "${sample_ids}_featureCounts.txt.summary", emit: read_counts_summary

    script:
    """
    if [ ${single_pair_end} -eq 2 ]; then
        featureCounts \
            -T $task.cpus \
            -p \
            --countReadPairs \
            -B \
            -C \
            -a $annotation_gtf_file \
            -o ${sample_ids}_featureCounts.txt \
            -t exon \
            -g gene_id \
            $aligned_bam
    else
        featureCounts \
            -T $task.cpus \
            -a $annotation_gtf_file \
            -o ${sample_ids}_featureCounts.txt \
            -t exon \
            -g gene_id \
            $aligned_bam
    fi
    """
}

// 22 Summarize statistics of assigned reads and counted
process _22_Read_count_mapped_stats {

    cpus params.threads
    tag "Counts_assigned_reads_summary" 
 
    publishDir "$params.run_name/06_gene_count", mode: "copy"

    input:
    path counts_summaries_stats

    output:
    path "Counts_stats_summary.csv", emit: counts_stats

    script:
    """
    echo "Run,Assigned_reads" > Counts_stats_summary.csv

    for summary in ${counts_summaries_stats}; do
        sample=\$(basename \$summary | sed 's/_featureCounts.txt.summary//')
        assigned=\$(awk '\$1 == "Assigned" {print \$2}' \$summary)
        echo "\$sample,\$assigned" >> Counts_stats_summary.csv
    done
    """
}

// 23 Merge statistic results into one table
process _23_Merge_stats {

    cpus params.threads
    tag "Summary_table"

    publishDir "$params.run_name", mode: "copy"

    input:
    val rows 

    output:
    path "Run_summary_table.csv"

    script:
    """
    echo "Run,Raw_reads,Avg_read_length,Library_name,Library_strategy,Platform,BioProject,Scientific_name,Sample_name,Sample_Control_manual_curated,Number_filtered_reads,Uniquely_mapped_perc,Multi_mapped_perc,Unmapped_short_perc,Unmapped_other_perc" > Run_summary_table.csv

    echo "$rows" >> "Run_summary_table.csv"
    """ 
}

// 24 Build one gene count matrix per BioProject
process _24_Gene_count_matrix_by_project {

    cpus params.threads
    tag "$bioproject-$plant_species"

    publishDir "$params.run_name/07_differential_expression/deseq2_gene_count_matrices", mode: "copy"

    input:
    tuple val(plant_species), val(bioproject), path(read_count_files)

    output:
    tuple val(plant_species), val(bioproject), path("${bioproject}_${plant_species}_gene_count_matrix.csv"), emit: gene_count_matrix

    script:
    """
    first=1

    for counts in ${read_count_files}; do
        sample=\$(basename \$counts | sed 's/_featureCounts.txt//')

        if [ \$first -eq 1 ]; then
            echo "Geneid,\$sample" > ${bioproject}_${plant_species}_gene_count_matrix.csv
            awk 'BEGIN {OFS=","} /^#/ {next} \$1 == "Geneid" {next} {print \$1,\$NF}' \$counts >> ${bioproject}_${plant_species}_gene_count_matrix.csv
            first=0
        else
            echo "\$sample" > \${sample}.counts.tmp
            awk 'BEGIN {OFS=","} /^#/ {next} \$1 == "Geneid" {next} {print \$NF}' \$counts >> \${sample}.counts.tmp
            paste -d',' ${bioproject}_${plant_species}_gene_count_matrix.csv \${sample}.counts.tmp > ${bioproject}_${plant_species}_gene_count_matrix.tmp
            mv ${bioproject}_${plant_species}_gene_count_matrix.tmp ${bioproject}_${plant_species}_gene_count_matrix.csv
            rm \${sample}.counts.tmp
        fi
    done
    """
}

// 25 Build one DESeq2 metadata table per BioProject and plant species
process _25_DESeq2_metadata_by_project {

    cpus params.threads
    tag "$bioproject-$plant_species"

    publishDir "$params.run_name/07_differential_expression/deseq2_metadata_tables", mode: "copy"

    input:
    tuple val(plant_species), val(bioproject), val(metadata_rows)

    output:
    tuple val(plant_species), val(bioproject), path("${bioproject}_${plant_species}_deseq2_metadata.csv"), emit: deseq2_metadata

    script:
    """
    echo "Run,Sample_Control_manual_curated" > ${bioproject}_${plant_species}_deseq2_metadata.csv
    echo "$metadata_rows" >> ${bioproject}_${plant_species}_deseq2_metadata.csv
    """
}

// 26 Run DESeq2 for each BioProject and plant species
process _26_Run_DESeq2 {

    cpus params.threads
    tag "$bioproject-$plant_species"

    publishDir "$params.run_name/07_differential_expression/deseq2_results", mode: "copy"

    input:
    tuple val(plant_species), val(bioproject), path(gene_count_matrix), path(deseq2_metadata)

    output:
    tuple val(plant_species), val(bioproject), path("${bioproject}_${plant_species}_deseq2_results.csv"), emit: deseq2_results
    path "${bioproject}_${plant_species}_deseq2_normalized_counts.csv", emit: deseq2_normalized_counts
    path "${bioproject}_${plant_species}_deseq2_status.txt", emit: deseq2_status

    script:
    """
    Rscript /data2/home/joaoc/Projects/VIGSt/VIGSt_nextflow/Scripts/run_deseq2.R \
        "${gene_count_matrix}" \
        "${deseq2_metadata}" \
        "${bioproject}_${plant_species}"
    """
}

// 27 Extract predicted protein sequences from plant genome annotations
process _27_Extract_plant_proteins {

    cpus params.threads
    tag "$plant_species"

    publishDir "$params.run_name/08_orthology/proteins", mode: "copy"

    input:
    tuple val(plant_species), path(genome_fasta_file), path(annotation_gtf_file)

    output:
    path "${plant_species}.cds.gtf"
    tuple val(plant_species), path("${plant_species}_gene_to_protein.tsv"), emit: gene_to_protein_map
    tuple val(plant_species), path("${plant_species}.faa"), emit: plant_proteins

    script:
    """
    awk -v map_out="${plant_species}_gene_to_protein.tsv" -v protein_prefix="${plant_species}|" 'BEGIN {FS=OFS="\\t"; print "Geneid\\tProteinID" > map_out}
        /^#/ {next}
        \$3 == "CDS" && (\$7 == "+" || \$7 == "-") {
            gene = ""; tx = ""; prot = ""

            if (match(\$9, /gene_id "[^"]+"/)) {
                gene = substr(\$9, RSTART + 9, RLENGTH - 10)
            }
            if (match(\$9, /transcript_id "[^"]+"/)) {
                tx = substr(\$9, RSTART + 15, RLENGTH - 16)
            }
            if (match(\$9, /protein_id "[^"]+"/)) {
                prot = substr(\$9, RSTART + 12, RLENGTH - 13)
            }

            if (tx != "" && tx !~ /^unassigned_transcript_/) {
                if (gene == "") gene = tx
                if (prot != "") tx = prot

                gene_out = gene
                protein_out = tx
                gsub(/\\./, "X", protein_out)

                gsub(/[^A-Za-z0-9_.:-]/, "_", gene)
                gsub(/[^A-Za-z0-9_.:-]/, "_", tx)

                print gene_out, protein_prefix protein_out >> map_out

                \$9 = "gene_id \\"" gene "\\"; transcript_id \\"" tx "\\";"
                print
            }
        }
    ' $annotation_gtf_file > ${plant_species}.cds.gtf

    gffread ${plant_species}.cds.gtf \
        -g $genome_fasta_file \
        -y ${plant_species}.raw.faa

    awk -v prefix="${plant_species}|" '
        /^>/ {
            sub(/^>/, ">" prefix)
            sub(/[ \\t].*/, "", \$0)
        }
        { print }
    ' ${plant_species}.raw.faa > ${plant_species}.aa.faa

    sed 's/\\./X/g' ${plant_species}.aa.faa > ${plant_species}.faa
    """
}

// 28 Infer orthogroups across plant species with OrthoFinder
process _28_Run_OrthoFinder {

   cpus params.threads
   tag "OrthoFinder"

   publishDir "$params.run_name/08_orthology", mode: "copy"

   input:
   path protein_fastas

   output:
   path "orthofinder_input"
   path "OrthoFinder_results", emit: orthofinder_results

   script:
   """
   mkdir -p orthofinder_input
   cp ${protein_fastas} orthofinder_input

   primary_transcript orthofinder_input

   orthofinder \
       -f orthofinder_input/primary_transcripts \
       -t $task.cpus \
       -a $task.cpus \
       -o OrthoFinder_results
   """
}

// 29 Add orthogroup IDs to DESeq2 results
process _29_Annotate_DESeq2_with_orthogroups {

    cpus params.threads
    tag "$bioproject-$plant_species"

    publishDir "$params.run_name/08_orthology/deseq2_orthogroups", mode: "copy"

    input:
    tuple val(plant_species), val(bioproject), path(deseq2_results), path(gene_to_protein_map), path(orthofinder_results)

    output:
    path "${bioproject}_${plant_species}_deseq2_orthogroups.csv", emit: deseq2_orthogroups

    script:
    """
    Rscript /data2/home/joaoc/Projects/VIGSt/VIGSt_nextflow/Scripts/annotate_deseq2_with_orthogroups.R \
        "${deseq2_results}" \
        "${gene_to_protein_map}" \
        "${orthofinder_results}" \
        "${plant_species}" \
        "${bioproject}_${plant_species}"
    """
}

// 30 Annotate plant proteins with eggNOG-mapper
process _30_Run_eggNOG_mapper {

    cpus params.threads
    tag "$plant_species"

    publishDir "$params.run_name/09_eggnog/eggnog_mapper", mode: "copy"

    input:
    tuple val(plant_species), path(protein_fasta)

    output:
    tuple val(plant_species), path("${plant_species}.emapper.annotations"), emit: eggnog_annotations

    script:
    """
    emapper.py \
        -i $protein_fasta \
        --itype proteins \
        --cpu $task.cpus \
        --output ${plant_species} \
        --output_dir . \
        --override
    """
}

// 31 Add eggNOG functional annotations to DESeq2 results
process _31_Annotate_DESeq2_with_eggNOG {

    cpus params.threads
    tag "$bioproject-$plant_species"

    publishDir "$params.run_name/09_eggnog/deseq2_eggnog", mode: "copy"

    input:
    tuple val(plant_species), val(bioproject), path(deseq2_results), path(gene_to_protein_map), path(eggnog_annotations)

    output:
    tuple val(plant_species), val(bioproject), path("${bioproject}_${plant_species}_deseq2_eggnog_annotations.csv"), emit: deseq2_eggnog_annotations

    script:
    """
    Rscript /data2/home/joaoc/Projects/VIGSt/VIGSt_nextflow/Scripts/annotate_deseq2_with_eggnog.R \
        "${deseq2_results}" \
        "${gene_to_protein_map}" \
        "${eggnog_annotations}" \
        "${plant_species}" \
        "${bioproject}_${plant_species}"
    """
}

// 32 Add OrthoFinder orthogroup IDs to per-study DESeq2 eggNOG annotations
process _32_Add_orthogroups_to_eggNOG_DESeq2 {

    cpus params.threads
    tag "$bioproject-$plant_species"

    publishDir "$params.run_name/09_eggnog/deseq2_eggnog_orthogroups", mode: "copy"

    input:
    tuple val(plant_species), val(bioproject), path(deseq2_eggnog_annotations), path(gene_to_protein_map), path(orthogroups_tsv)

    output:
    path "${bioproject}_${plant_species}_deseq2_eggnog_orthogroups.csv", emit: deseq2_eggnog_orthogroups

    script:
    """
    Rscript /data2/home/joaoc/Projects/VIGSt/VIGSt_nextflow/Scripts/add_orthogroups_to_eggnog_annotations.R \
        "${deseq2_eggnog_annotations}" \
        "${gene_to_protein_map}" \
        "${orthogroups_tsv}" \
        "${plant_species}" \
        "${bioproject}_${plant_species}"
    """
}

// 33 Summarize DEG log2 fold changes by orthogroup across all studies
process _33_Create_orthogroup_log2fc_summary {

    cpus params.threads
    tag "Orthogroup_log2FoldChange_summary"

    publishDir "$params.run_name/10_orthogroup_log2fc_summary", mode: "copy"

    input:
    path deseq2_eggnog_orthogroup_files

    output:
    path "orthogroup_log2FoldChange_summary.csv", emit: orthogroup_log2fc_summary

    script:
    """
    Rscript /data2/home/joaoc/Projects/VIGSt/VIGSt_nextflow/Scripts/create_orthogroup_log2fc_summary.R \\
        orthogroup_log2FoldChange_summary.csv \\
        ${deseq2_eggnog_orthogroup_files}
    """
}

// 34 Identify upregulated and downregulated orthogroups across plant species
process _34_Find_cross_plant_gene_expression_patterns {

    cpus params.threads
    tag "Cross_plant_gene_expression_patterns"

    publishDir "$params.run_name/10_orthogroup_log2fc_summary", mode: "copy"

    input:
    path orthogroup_log2fc_summary

    output:
    path "orthogroup_cross_plant_gene_expression_patterns.csv", emit: cross_plant_regulation_patterns

    script:
    """
    Rscript /data2/home/joaoc/Projects/VIGSt/VIGSt_nextflow/Scripts/find_cross_plant_gene_expression_patterns.R \\
        "${orthogroup_log2fc_summary}" \\
        orthogroup_cross_plant_gene_expression_patterns.csv
    """
}

// 35 Run TopGO enrichment for cross-plant expression-pattern orthogroups
process _35_Enrich_cross_plant_orthogroups {

    cpus params.threads
    tag "Cross_plant_orthogroup_TopGo_enrichment"

    publishDir "$params.run_name/11_orthogroup_functional_enrichment", mode: "copy"

    input:
    path cross_plant_patterns
    path deseq2_eggnog_orthogroup_files
    path topgo_script

    output:
    path "topgo_go_enrichment.csv", emit: topgo_enrichment
    path "cross_plant_patterns_with_topgo.csv", emit: patterns_with_topgo

    script:
    """
    Rscript ${topgo_script} \\
        "${cross_plant_patterns}" \\
        topgo_go_enrichment.csv \\
        cross_plant_patterns_with_topgo.csv \\
        ${deseq2_eggnog_orthogroup_files}
    """
}

// Workflow block
workflow {
    _01_Query_SRA_IDs(params.query)

    _04_Download_sample_project_metadata(params.query)
    
    _03_Exclude_studies(_01_Query_SRA_IDs.out.run_info, channel.fromPath("Non_NCBI_Data_TranscriptomicSamples_PlantGenomes/excluded_studies.csv"))

    _02_Exclude_samples(_03_Exclude_studies.out.run_info_studies_clean, channel.fromPath("Non_NCBI_Data_TranscriptomicSamples_PlantGenomes/excluded_samples.csv"))

    Sample_meta_ch = _02_Exclude_samples.out.run_info_simplified_clean 
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

    _05_Download_raw_fastq_files(Sample_meta_ch)  
    
    Short_fastqc_ch = _05_Download_raw_fastq_files.out.raw_fastq_files
        .filter {fastq_short_file -> fastq_short_file.toString().contains("/Short_reads/")}
        .flatten()

        .map {fastq_short ->
            def bioproject = fastq_short.parent.name

            tuple(bioproject, fastq_short)}

    _06_FastQC_raw_samples(Short_fastqc_ch)

    Multiqc_ch = _06_FastQC_raw_samples.out.fastqc_zips
        .map {fastqc_zip_file -> 
            def bioproject = fastqc_zip_file.parent.name

            tuple(bioproject, fastqc_zip_file)}
            .groupTuple()

    _07_MultiQC_raw_samples(Multiqc_ch)

    Long_nanoplot_ch = _05_Download_raw_fastq_files.out.raw_fastq_files
        .filter {fastq_long_file -> fastq_long_file.toString().contains("/Long_reads/")}
        .flatten()

        .map {fastq_long ->
            def bioproject = fastq_long.parent.name
            def sample_ids = fastq_long.simpleName

            tuple(bioproject, sample_ids, fastq_long)}

    _08_NanoPlot_raw_samples(Long_nanoplot_ch)

     Filtered_short_fastq_ch = Short_fastqc_ch
        .map {_bioproject, fastq_short ->
            def bioproject = fastq_short.parent.name
            def sample_ids = fastq_short.simpleName.replaceAll(/(_[12])$/, '')
            def fastq_short_reads = fastq_short

            tuple(bioproject, sample_ids, fastq_short_reads)}
        .groupTuple(by:[0,1])

    _09_Fastp_filtering(Filtered_short_fastq_ch)

    Filtered_long_fastq_ch = Long_nanoplot_ch
    
    _10_Fastplong_filtering(Filtered_long_fastq_ch)

    Viroid_refseq_ch = _11_Download_viroid_genomes()

    Viroid_index_ch = _12_Build_viroid_index(Viroid_refseq_ch.viroid_fasta)
        
    _13_Map_to_viroid_genomes(_09_Fastp_filtering.out.filtered_short_fastq_files, Viroid_index_ch)

    _14_Map_to_viroid_genomes_stats(_09_Fastp_filtering.out.filtered_short_fastq_files, _13_Map_to_viroid_genomes.out.aligned_bam)

    _15_Merge_viroid_csvs(_14_Map_to_viroid_genomes_stats.out.virod_stats.collect())

    Plant_genomes_refseq_ch = channel.of(
        tuple("Solanum_tuberosum", "GCF/000/226/075/GCF_000226075.1_SolTub_3.0"),
        tuple("Solanum_lycopersicum", "GCF/036/512/215/GCF_036512215.1_SLM_r2.1"),
        tuple("Cucumis_sativus", "GCF/000/004/075/GCF_000004075.3_Cucumber_9930_V3"),
        tuple("Humulus_lupulus", "GCF/963/169/125/GCF_963169125.1_drHumLupu1.1"))
        //tuple("Citrus_x_aurantiifolia","GCA/052/148/995/GCA_052148995.1_Mxlime_USDA_v1"),
        //tuple("Citrus", "GCF/022/201/045/GCF_022201045.2_DVS_A1.0"),
        //tuple("Prunus_persica", "GCF/000/346/465/GCF_000346465.2_Prunus_persica_NCBIv2"),
        //tuple("Vitis_vinifera", "GCF/030/704/535/GCF_030704535.1_ASM3070453v1"),
        //tuple("Cocos_nucifera", "GCA/008/124/465/GCA_008124465.1_ASM812446v1"),
        //tuple("Malus_domestica", "GCF/042/453/785/GCF_042453785.1_GDT2T_hap1"))

    .map  {plant_species, assembly ->
        def base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/${assembly}"
        def assembly_name = assembly.tokenize("/").last()
        def genome_url = "${base_url}/${assembly_name}_genomic.fna.gz"
        def gtf_url = "${base_url}/${assembly_name}_genomic.gtf.gz"

        tuple(plant_species, genome_url, gtf_url)}

    _16_Download_plant_genomes(Plant_genomes_refseq_ch)

    Plant_genomes_downloaded_ch = channel.of()
        //tuple("Citrus_medica", file("Non_NCBI_Data_TranscriptomicSamples_PlantGenomes/genomes/Citrus_medica/genome.fna"), file("Non_NCBI_Data_TranscriptomicSamples_PlantGenomes/genomes/Citrus_medica/annotation.gff3")),
        //tuple("Rubus_idaeus", file("Non_NCBI_Data_TranscriptomicSamples_PlantGenomes/genomes/Rubus_idaeus/genome.fna"), file("Non_NCBI_Data_TranscriptomicSamples_PlantGenomes/genomes/Rubus_idaeus/annotation.gff3")),
        //tuple("Chrysanthemum_x_morifolium", file("Non_NCBI_Data_TranscriptomicSamples_PlantGenomes/genomes/Chrysanthemum_x_morifolium/genome.fna"), file("Non_NCBI_Data_TranscriptomicSamples_PlantGenomes/genomes/Chrysanthemum_x_morifolium/annotation.gff3")))

    _17_Downloaded_plant_genomes(Plant_genomes_downloaded_ch)

    Plant_genomes_ch = _16_Download_plant_genomes.out.mix(_17_Downloaded_plant_genomes.out)

    _18_Build_plant_index(Plant_genomes_ch)

    Genomes_metadata_ch = _02_Exclude_samples.out.run_info_simplified_clean
        .splitCsv(header:true)
        .map {row -> 

        def sample_ids = row.Run.trim()
        def plant_species = row.ScientificName.trim().replaceAll(' ', '_')    
        tuple(sample_ids, plant_species)}

    Samples_ch = _09_Fastp_filtering.out.filtered_short_fastq_files
        .map {bioproject, sample_id, reads ->
        tuple(sample_id, bioproject, reads)}

    Samples_with_species_ch = Samples_ch.join(Genomes_metadata_ch)
        .map {bioproject, sample_id, reads, plant_species ->
        tuple(plant_species, sample_id, bioproject, reads)}

    Plant_index_ch = _18_Build_plant_index.out.plant_index

    Samples_with_index_ch = Samples_with_species_ch
        .map {plant_species, sample_id, bioproject, reads ->
        tuple(plant_species, sample_id, bioproject, reads, reads.size())}
        .combine(Plant_index_ch, by: 0)
    
    _19_Map_to_plant_genomes(Samples_with_index_ch)

    _20_Map_to_plant_genomes_stats(_19_Map_to_plant_genomes.out.log_map_stats.collect())

    Run_info_simplified_ch = _02_Exclude_samples.out.run_info_simplified_clean
        .splitCsv(header:true)
        .map {row -> tuple(row.Run, row.spots, row.avgLength, row.LibraryName, row.LibraryStrategy, row.Platform, row.BioProject, row.ScientificName, row.SampleName)}

    Sample_control_curated_ch = channel.fromPath("Non_NCBI_Data_TranscriptomicSamples_PlantGenomes/viroid_RNA_Seq_runs_info_sample_control_curated.csv")
        .splitCsv(header:true)
        .map {row -> tuple(row.Run, row.Sample_Control_manual_curated)}

    //viroid_map_stats_ch = _15_Merge_viroid_csvs

    Plant_map_stats_ch =_20_Map_to_plant_genomes_stats.out.plant_map_stats
        .splitCsv(header:true)
        .map {row -> tuple(row.Run, row.Number_filtered_reads, row.Uniquely_mapped_perc, row.Multi_mapped_perc, row.Unmapped_short_perc, row.Unmapped_other_perc)}

    Plant_only_annotation_ch = Plant_genomes_ch
        .map {plant_species, _genome, gtf ->
        tuple(plant_species, gtf)}

    Read_count_ch =_19_Map_to_plant_genomes.out.aligned_bam
        .combine(Plant_only_annotation_ch, by: 0)

    _21_Read_count_mapped(Read_count_ch)

    _22_Read_count_mapped_stats(_21_Read_count_mapped.out.read_counts_summary.collect())

    Join_summary_tables_ch = Run_info_simplified_ch.join(Sample_control_curated_ch, by: 0, remainder:true).join(Plant_map_stats_ch, by: 0, remainder:true)
        .map {row -> row.join(",")}
        .collect()
        .map {rows -> rows.join('\n')}
    
    _23_Merge_stats(Join_summary_tables_ch)

    Read_counts_by_project_ch = _21_Read_count_mapped.out.read_counts
        .map {plant_species, bioproject, _sample_ids, read_count_file ->
        tuple(plant_species, bioproject, read_count_file)}
        .groupTuple(by: [0,1])

    _24_Gene_count_matrix_by_project(Read_counts_by_project_ch)

    DESeq2_metadata_by_project_ch = _21_Read_count_mapped.out.read_counts
        .map {plant_species, bioproject, sample_ids, _read_count_file ->
        tuple(sample_ids, plant_species, bioproject)}
        .join(Sample_control_curated_ch)
        .map {sample_ids, plant_species, bioproject, sample_control ->
        tuple(plant_species, bioproject, "${sample_ids},${sample_control}")}
        .groupTuple(by: [0,1])
        .map {plant_species, bioproject, metadata_rows ->
        tuple(plant_species, bioproject, metadata_rows.join('\n'))}

    _25_DESeq2_metadata_by_project(DESeq2_metadata_by_project_ch)

    DESeq2_input_ch = _24_Gene_count_matrix_by_project.out.gene_count_matrix
        .join(_25_DESeq2_metadata_by_project.out.deseq2_metadata, by: [0,1])

    _26_Run_DESeq2(DESeq2_input_ch)

    _27_Extract_plant_proteins(Plant_genomes_ch)

    _28_Run_OrthoFinder(_27_Extract_plant_proteins.out.plant_proteins
        .map{_plant_species, protein_fasta -> protein_fasta}
        .collect())

    _30_Run_eggNOG_mapper(_27_Extract_plant_proteins.out.plant_proteins)

    DESeq2_eggnog_input_ch = _26_Run_DESeq2.out.deseq2_results
        .combine(_27_Extract_plant_proteins.out.gene_to_protein_map, by: 0)
        .combine(_30_Run_eggNOG_mapper.out.eggnog_annotations, by: 0)

    _31_Annotate_DESeq2_with_eggNOG(DESeq2_eggnog_input_ch)

    DESeq2_eggnog_orthogroup_input_ch = _31_Annotate_DESeq2_with_eggNOG.out.deseq2_eggnog_annotations
        .combine(_27_Extract_plant_proteins.out.gene_to_protein_map, by: 0)
        .combine(channel.fromPath("Viroid_Viridiplantae_transcriptomic_02/08_orthology/OrthoFinder_results/Results_Jul28/Orthogroups/Orthogroups.tsv"))

    _32_Add_orthogroups_to_eggNOG_DESeq2(DESeq2_eggnog_orthogroup_input_ch)

    _33_Create_orthogroup_log2fc_summary(
        _32_Add_orthogroups_to_eggNOG_DESeq2.out.deseq2_eggnog_orthogroups.collect())

    Orthogroup_pattern_input_ch = channel.fromPath("Viroid_Viridiplantae_transcriptomic_02/09_eggnog/deseq2_eggnog_orthogroups_run/orthogroup_log2FoldChange_summary.csv")
    _34_Find_cross_plant_gene_expression_patterns(Orthogroup_pattern_input_ch)

    DESeq2_eggnog_orthogroups_ch = channel.fromPath("Viroid_Viridiplantae_transcriptomic_02/09_eggnog/deseq2_eggnog_orthogroups_run/*_deseq2_eggnog_orthogroups.csv").collect()
    TopGO_script_ch = channel.fromPath("Scripts/run_topgo_cross_plant_orthogroups.R")
    _35_Enrich_cross_plant_orthogroups(_34_Find_cross_plant_gene_expression_patterns.out.cross_plant_regulation_patterns, DESeq2_eggnog_orthogroups_ch, TopGO_script_ch)

}
