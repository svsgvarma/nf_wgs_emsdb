#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */

// Primary input (file of input files, one per line)
params.reads_bam = "${params.input_dir}/data/sample_bams.txt"

// Accessory files
params.reference        = "${params.input_dir}/data/ref/ref.fasta"
params.reference_index  = "${params.input_dir}/data/ref/ref.fasta.fai"
params.reference_dict   = "${params.input_dir}/data/ref/ref.dict"
params.intervals        = "${params.input_dir}/data/ref/intervals.bed"

// Base name for final output file
params.cohort_name = "family_trio"

include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk/haplotypecaller/main.nf'
include { GATK_JOINTGENOTYPING } from './modules/gatk/jointgenotyping/main.nf'
include { MULTIQC } from './modules/multiqc/main.nf'

workflow {
    // Create input channel from a text file listing input file paths
    reads_ch = Channel.fromPath(params.reads_bam).splitText()

    // Load the file paths for the accessory files (reference and intervals)
    ref_file        = file(params.reference)
    ref_index_file  = file(params.reference_index)
    ref_dict_file   = file(params.reference_dict)
    intervals_file  = file(params.intervals)

    // Create index file for input BAM file
    SAMTOOLS_INDEX(reads_ch)

    // Call variants from the indexed BAM file
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
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
