process PREPARE_VCF {
    tag "${cohort}"
    publishDir "${params.outdir}/prepared", mode: 'copy', overwrite: true
    cpus 4
    memory '8 GB'
    input:
        path vcf
        path fasta
        val cohort
        val keep_filtered
    output:
        tuple path("${cohort}.prepared.vcf.gz"), path("${cohort}.prepared.vcf.gz.tbi"), emit: vcf
        path "${cohort}.prepare_*", emit: stats
    script:
    """
    prepare_vcf.sh ${vcf} ${fasta} ${cohort} ${params.normalize} ${params.hard_filter} ${keep_filtered} ${task.cpus}
    """
}