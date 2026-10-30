process SNP_QC_REPORT {
    publishDir "${params.outdir}/qc_report", mode: 'copy', overwrite: true
    cpus params.qc_report_cpus
    memory '16 GB'
    input:
        path qc_files
        path joint_vcf
        val cohort_name
    output:
        path "${cohort_name}.snp_calling_qc_report.html", emit: html
        path "${cohort_name}.*.tsv", emit: tsv
    script:
    def chrom_pattern = params.qc_chrom_pattern ? "--chrom-pattern '${params.qc_chrom_pattern}'" : ''
    def regions = params.intervals ? file(params.intervals).name :
        ('all contigs' + (params.contig_regex ? " matching ${params.contig_regex}" : '') +
         (params.min_contig_length ? " >= ${params.min_contig_length} bp" : ''))
    """
    case "${joint_vcf}" in
        *.gz) zcat ${joint_vcf} > joint_for_report.vcf; VCF=joint_for_report.vcf ;;
        *) VCF=${joint_vcf} ;;
    esac
    snp_qc_report.py --qc-dir . --vcf \$VCF --cohort ${cohort_name} --outdir . \\
        --reference '${file(params.reference).name}' --intervals '${regions}' \\
        --min-mapping-rate ${params.qc_min_mapping_rate} --min-mean-depth ${params.qc_min_mean_depth} \\
        --max-duplicate-rate ${params.qc_max_duplicate_rate} --max-missing-rate ${params.qc_max_missing_rate} \\
        --max-contigs ${params.qc_max_contigs} --many-samples ${params.qc_many_samples} \\
        --threads ${task.cpus} ${chrom_pattern}
    rm -f joint_for_report.vcf
    """
}