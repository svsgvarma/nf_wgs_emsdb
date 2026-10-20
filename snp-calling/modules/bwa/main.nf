/*
 * Build BWA index with custom prefix, keeping reference genome name
 */
process BWA_INDEX {
    // Publish index files alongside the reference FASTA in the input S3 bucket
    // so they are available for future runs without rebuilding
    publishDir "${file(params.reference).parent}", mode: 'copy'

    input:
        path ref_fasta

    output:
        tuple val("bwa_index"), path("${ref_fasta}.*")

    script:
    // Default bwa index (no -p) produces ref.fasta.amb, .ann, .bwt, .pac, .sa
    """
    bwa index "$ref_fasta"

    """
}

/*
 * Align paired-end reads using BWA MEM with read group and threads
 */

params.bwa_threads = 8

process BWA_ALIGN {
    // SAM is never written to disk — BWA output is piped directly into samtools view to produce BAM.
    // This avoids staging a 100-200 GB SAM file from S3 in a downstream SAMTOOLS_CONVERT step,
    // which previously caused [Errno 12] Cannot allocate memory on the Batch instance.
    // SAMTOOLS_CONVERT is therefore no longer needed and has been removed from the workflow.
    // No need to publish to S3; disabling publishDir saves storage and upload cost.
    // publishDir params.outdir, mode: 'copy'
    cpus params.bwa_threads ?: 8
    memory '16 GB'

    input:
        tuple val(sample_id), path(reads1), path(reads2)
        // val index_prefix does not stage files from S3 — Nextflow only copies the string value, not the files.
        // Use path so all BWA index files (.amb .ann .bwt .pac .sa) are staged to the work directory.
        path bwa_index_files

    output:
        tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    index_prefix=\$(ls *.amb | sed 's/\\.amb//')
    bwa mem -t ${task.cpus} -R "@RG\\tID:${sample_id}\\tSM:${sample_id}" \$index_prefix $reads1 $reads2 | samtools view -bS - > ${sample_id}.bam
    """
}

