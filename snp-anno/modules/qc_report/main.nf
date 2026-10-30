process ANNOTATION_QC_REPORT {
    tag "${cohort}"
    publishDir "${params.outdir}/qc_report", mode: 'copy', overwrite: true
    cpus 1
    memory '8 GB'
    input:
        tuple path(vcf), path(vcf_index)
        path stats_files
        val cohort
        val genome
        val input_vcf
    output:
        path "${cohort}.annotation_qc_report.html", emit: html
        path "${cohort}.*.tsv", emit: tsv
    script:
    """
    annotation_qc_report.py --vcf ${vcf} --stats-dir . --cohort ${cohort} --genome ${genome} \\
        --input-vcf '${input_vcf}' --snpeff-summary '../annotation/${cohort}.snpeff_summary.html' \\
        --max-contigs ${params.qc_max_contigs} --many-samples ${params.qc_many_samples} \\
        --top-genes ${params.qc_top_genes} --window ${(params.qc_window_mb * 1000000) as long} --outdir .
    """
}