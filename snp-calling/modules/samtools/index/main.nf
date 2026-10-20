/*
 * Convert SAM to BAM and index BAM files using samtools
 */

// SAMTOOLS_CONVERT has been removed from the workflow.
// BWA_ALIGN now pipes its output directly into samtools view, producing BAM without writing SAM to disk.
// A WGS SAM file is 100-200 GB per sample — staging it from S3 to the Batch instance caused
// [Errno 12] Cannot allocate memory. The pipe approach eliminates this entirely.
// This process is kept here for reference but is no longer called.
/*
 process SAMTOOLS_CONVERT {
    // Unsorted BAM is an intermediate file — replaced by the sorted BAM from SAMTOOLS_SORT.
    // No need to publish to S3; disabling publishDir saves storage and upload cost.
    // publishDir params.outdir, mode: 'copy'

    cpus 2
    memory '64 GB'

    input:
        tuple val(sample_id), path(sam_file)

    output:
        tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    samtools view -bS $sam_file > ${sample_id}.bam
    """
}
*/

process SAMTOOLS_SORT {
    // mode: 'symlink' does not work with S3 — symlinks are local filesystem constructs and cannot be created in S3. Using copy instead
    // publishDir params.outdir, mode: 'symlink'
    publishDir params.outdir, mode: 'copy'

    memory '16 GB'

    input:
        tuple val(sample_id), path(bam_file)

    output:
        tuple val(sample_id), path("${sample_id}.sorted.bam")

    script:
    """
    samtools sort -o ${sample_id}.sorted.bam $bam_file
    """
}

process SAMTOOLS_ADDRG {
    // mode: 'symlink' does not work with S3 — symlinks are local filesystem constructs and cannot be created in S3. Using copy instead
    // publishDir params.outdir, mode: 'symlink'
    publishDir params.outdir, mode: 'copy'

    memory '16 GB'

    input:
        tuple val(sample_id), path(sorted_bam)

    output:
        tuple val(sample_id), path("${sample_id}.rg.bam")

    script:
    """
    samtools addreplacerg \
        -r 'ID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA' \
        -o ${sample_id}.rg.bam \
        $sorted_bam
    """
}

process SAMTOOLS_INDEX {
    // mode: 'symlink' does not work with S3 — symlinks are local filesystem constructs and cannot be created in S3. Using copy instead
    // publishDir params.outdir, mode: 'symlink'
    publishDir params.outdir, mode: 'copy'

    memory '16 GB'

    input:
        tuple val(sample_id), path(bam_file)

    output:
        tuple val(sample_id), path(bam_file), path("${bam_file}.bai")

    script:
    """
    samtools index $bam_file
    """
}
