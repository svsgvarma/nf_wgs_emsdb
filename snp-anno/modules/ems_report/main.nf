process EMS_REPORT {
    tag "${cohort}"
    publishDir "${params.outdir}/ems_report", mode: 'copy', overwrite: true
    cpus params.ems_report_cpus
    memory '16 GB'
    input:
        tuple path(vcf), path(vcf_index)
        path alignment_qc, stageAs: 'opt_alignment_qc/*'
        path candidate_genes, stageAs: 'opt_candidate_genes/*'
        path homeolog_pairs, stageAs: 'opt_homeolog_pairs/*'
        val cohort
    output:
        path "${cohort}.ems_report.html", emit: html
        path "${cohort}.ems_*.tsv*", emit: tables
    script:
    def alignment_qc_arg = alignment_qc.toString().endsWith('NO_FILE') ? "''" : "'${alignment_qc}'"
    def candidate_genes_arg = candidate_genes.toString().endsWith('NO_FILE') ? "''" : "'${candidate_genes}'"
    def homeolog_pairs_arg = homeolog_pairs.toString().endsWith('NO_FILE') ? "''" : "'${homeolog_pairs}'"
    """
    set -o pipefail
    bcftools view -h ${vcf} > header.txt
    bcftools query -l ${vcf} > samples.txt
    MQ_FIELD=\$(grep -q '^##INFO=<ID=MQ,' header.txt && echo '%INFO/MQ' || echo '.')
    bcftools query -i 'TYPE="snp"' \\
        -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t\${MQ_FIELD}\\t%FILTER\\t%INFO/ANN[\\t%AD:%DP]\\n" \\
        ${vcf} > sites.tsv
    ems_report.py --sites sites.tsv --samples samples.txt --header header.txt --cohort ${cohort} --outdir . \\
        --wt-samples '${params.ems_wt_samples}' --mutant-samples '${params.ems_mutant_samples}' \\
        --pool-sizes '${params.ems_pool_sizes}' --families '${params.ems_families}' \\
        --min-alt ${params.ems_min_alt} --min-dp ${params.ems_min_dp} --min-af ${params.ems_min_af} \\
        --max-dp ${params.ems_max_dp} --max-dp-factor ${params.ems_max_dp_factor} \\
        --min-mq ${params.ems_min_mq} --min-qual ${params.ems_min_qual} --wt-min-alt ${params.ems_wt_min_alt} \\
        --hom-af ${params.ems_hom_af} --ploidy ${params.ems_ploidy} --min-purity ${params.ems_min_purity} \\
        --canonical-only ${params.ems_canonical_only} --alignment-qc ${alignment_qc_arg} \\
        --candidate-genes ${candidate_genes_arg} --homeolog-pairs ${homeolog_pairs_arg} \\
        --window ${(params.qc_window_mb * 1000000) as long} --max-contigs ${params.qc_max_contigs} \\
        --many-samples ${params.qc_many_samples} --top-genes ${params.qc_top_genes} --threads ${task.cpus}
    rm -f sites.tsv
    """
}