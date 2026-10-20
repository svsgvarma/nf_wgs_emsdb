#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */

// Primary input (file of input files, one per line)

params.reads_fastq = "${params.input_dir}/data/fastq/NIST_Illumina_2x250bps_*/*_R{1,2}_*.fastq.gz"
//params.reads_bam = "${params.input_dir}/data/sample_bams.txt"

// Output directory commented out because this is set in the nextflow.config file when running in AWS
//params.outdir = "results" 

// Accessory files
params.reference        = "${params.input_dir}/data/ref/GRCh38_GIAB/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
params.reference_index  = "${params.input_dir}/data/ref/GRCh38_GIAB/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai"
params.reference_index_bwa  = "${params.input_dir}/data/ref/GRCh38_GIAB/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
params.reference_dict   = "${params.input_dir}/data/ref/GRCh38_GIAB/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict"
params.intervals        = "${params.input_dir}/data/ref/GRCh38_GIAB/GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2.bed"

// Base name for final output file
params.cohort_name = "family_trio"

include { BWA_INDEX; BWA_ALIGN } from './modules/bwa/main.nf'
// SAMTOOLS_CONVERT removed — BWA_ALIGN now pipes directly to BAM via samtools view
include { SAMTOOLS_SORT; SAMTOOLS_ADDRG; SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk/haplotypecaller/main.nf'
include { GATK_JOINTGENOTYPING } from './modules/gatk/jointgenotyping/main.nf'
include { MULTIQC } from './modules/multiqc/main.nf'

workflow {
    // Create input channel from a glob pattern for paired-end reads
    reads_ch = Channel.fromFilePairs(params.reads_fastq, flat: true)

    // Load reference FASTA
    ref_file = file(params.reference)
    ref_index_file = file(params.reference_index)   
    ref_dict_file  = file(params.reference_dict)
    intervals_file = file(params.intervals)
    
    // Load pre-built BWA index files from S3 (must exist alongside the reference FASTA)
    // To build: bwa index <ref.fasta>  then upload the 5 files (.amb .ann .bwt .pac .sa) to S3
    bwa_ext = ['amb', 'ann', 'bwt', 'pac', 'sa']
    bwa_ext.each { ext ->
        if (!file("${params.reference_index_bwa}.${ext}").exists())
            error "Missing BWA index file: ${params.reference_index_bwa}.${ext} — run 'bwa index' and upload to S3"
    }
    bwa_index_files = Channel.fromPath(
        bwa_ext.collect { ext -> "${params.reference_index_bwa}.${ext}" }
    ).collect()

    // Run Bwa alignment — output is BAM (SAM is never written to disk, piped internally via samtools view)
    bwa_out = BWA_ALIGN(reads_ch, bwa_index_files)

    // Sort BAM files
    sorted_bam_files_ch = SAMTOOLS_SORT(bwa_out)

    // Add read group to sorted BAM files
    //rg_bam_files_ch = SAMTOOLS_ADDRG(sorted_bam_files_ch)

    // Index sorted BAM files
    indexed_bam_ch = SAMTOOLS_INDEX(sorted_bam_files_ch)

    // Prepare input for GATK_HAPLOTYPECALLER: tuple of BAM and BAI
    gatk_input_ch = indexed_bam_ch.map { tuple -> [tuple[1], tuple[2]] }

    // Call variants from the indexed BAM file
    GATK_HAPLOTYPECALLER(
        gatk_input_ch,
        ref_file,
        ref_index_file,
        ref_dict_file,
        intervals_file
    )

    // Collect variant calling outputs across samples
    all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
    all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()

    // Combine GVCFs into a GenomicsDB data store and apply joint genotyping
    GATK_JOINTGENOTYPING(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name,
        ref_file,
        ref_index_file,
        ref_dict_file
    )

    // // Collect outputs for MultiQC
    // multiqc_input = Channel.merge(
    //     SAMTOOLS_INDEX.out,
    //     GATK_HAPLOTYPECALLER.out.vcf,
    //     GATK_JOINTGENOTYPING.out.vcf
    // )

    // // Run MultiQC to generate reports
    // MULTIQC(multiqc_input, params.outdir)
}

//##############  Cleanup step  ##############
// Cleanup step to remove log files, ./work/ sub-directories, and .nextflow/ folder after workflow completion
//############################################
 
workflow.onComplete {
    if (workflow.success) {
        println "Workflow succeeded. Running cleanup commands."

        // Execute cleanup commands
        def cleanupScript = """
        echo "Deleting work/ folder, .nextflow/ folder, and all log files."

        # Remove all directories and files inside the work/ folder
        if [ -d "work" ]; then
            rm -rf work/*
        fi

        # Remove the .nextflow/ folder
        if [ -d ".nextflow" ]; then
            rm -rf .nextflow/
        fi

        # Remove all log files in the current directory, including rotated ones
        find . -name "*.log*" -type f -delete

        # Remove all index log files in the /results directory
        #find . -name "*.bt2*" -type f -delete
        """
        
        // Execute the cleanup commands using bash
        def process = ["bash", "-c", cleanupScript].execute()
        process.waitFor()

        // Check for errors
        if (process.exitValue() != 0) {
            println "Cleanup script failed with error: ${process.err.text}"
        } else {
            println "Cleanup completed successfully."
        }
    } else {
        println "Workflow failed. Preserving work/ folder, .nextflow/ folder, and log files for debugging."
    }
}


/*
process ALIGN {
    cpus 8
    memory '32 GB'
    time '4h'

    input:
    // define your input here

    output:
    // define your output here

    script:
    """
    # your alignment command here
    """
}
*/
