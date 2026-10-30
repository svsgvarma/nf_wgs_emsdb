process SAMTOOLS_MARKDUP {
    tag "${sample_id}"
    publishDir params.outdir, mode: 'copy', pattern: '*.markdup.bam', enabled: params.publish_bams == true
    cpus 8
    memory { 16.GB * task.attempt }
    input:
        tuple val(sample_id), path(bam_file)
    output:
        tuple val(sample_id), path("${sample_id}.markdup.bam"), emit: bam
        path "${sample_id}.qc.markdup.txt", emit: metrics
    script:
    """
    set -o pipefail
    samtools fixmate -@ ${task.cpus} -m ${bam_file} - \\
        | samtools sort -@ ${task.cpus} -m 1G -T sort_tmp -o - - \\
        | samtools markdup -@ ${task.cpus} -f ${sample_id}.qc.markdup.txt - ${sample_id}.markdup.bam
    """
}