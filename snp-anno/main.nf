#!/usr/bin/env nextflow

params.vcf = null
params.cohort_name = "family_trio"
params.genome_fasta = null
params.genome_gff = null
params.genome_name = params.input_subdir
params.snpeff_codon_tables = ''
params.normalize = true
params.hard_filter = true
params.keep_filtered = false
params.snpeff_args = ''
params.snpsift_fields = 'CHROM POS REF ALT FILTER AF AC ANN[0].EFFECT ANN[0].IMPACT ANN[0].GENE ANN[0].GENEID ANN[0].FEATUREID ANN[0].BIOTYPE ANN[0].RANK ANN[0].HGVS_C ANN[0].HGVS_P ANN[0].ERRORS'
params.qc_max_contigs = 40
params.qc_many_samples = 24
params.qc_top_genes = 50
params.qc_window_mb = 1
params.ems = false
params.ems_wt_samples = ''
params.ems_mutant_samples = ''
params.ems_pool_sizes = ''
params.ems_families = ''
params.ems_min_alt = 3
params.ems_min_dp = 6
params.ems_min_af = 0
params.ems_max_dp = 0
params.ems_max_dp_factor = 3
params.ems_min_mq = 40
params.ems_min_qual = 0
params.ems_wt_min_alt = 2
params.ems_hom_af = 0.8
params.ems_ploidy = 2
params.ems_min_purity = 80
params.ems_canonical_only = true
params.ems_candidate_genes = ''
params.ems_homeolog_pairs = ''
params.ems_alignment_qc = "${params.snp_calling_results}/qc_report/${params.cohort_name}.sample_alignment_qc.tsv"
params.ems_report_cpus = 8

include { SNPEFF_BUILD } from './modules/snpeff/build/main.nf'
include { PREPARE_VCF } from './modules/bcftools/prepare/main.nf'
include { SNPEFF_ANNOTATE } from './modules/snpeff/annotate/main.nf'
include { SNPSIFT_EXTRACT } from './modules/snpsift/extract/main.nf'
include { ANNOTATION_QC_REPORT } from './modules/qc_report/main.nf'
include { EMS_REPORT } from './modules/ems_report/main.nf'

def optional_file(path) {
    return (path && file(path).exists()) ? file(path) : file("${projectDir}/assets/NO_FILE")
}

workflow {
    if (!params.genome_fasta || !params.genome_gff)
        error "Provide --genome_fasta and --genome_gff (the reference used for SNP calling and its GFF3/GTF annotation)"
    if (!params.genome_name)
        error "Provide --genome_name (defaults to input_subdir)"
    if (params.ems && !params.ems_wt_samples)
        error "--ems true needs --ems_wt_samples (comma-separated wild-type libraries used as background mask)"

    def default_vcf = "${params.snp_calling_results}/${params.cohort_name}.joint.vcf.gz"
    if (!params.vcf && !file(default_vcf).exists() && file(default_vcf - '.gz').exists())
        default_vcf = default_vcf - '.gz'
    vcf_path = params.vcf ?: default_vcf
    vcf = file(vcf_path, checkIfExists: true)
    fasta = file(params.genome_fasta, checkIfExists: true)
    gff = file(params.genome_gff, checkIfExists: true)
    keep_filtered = params.keep_filtered || params.ems

    SNPEFF_BUILD(fasta, gff, params.genome_name)
    PREPARE_VCF(vcf, fasta, params.cohort_name, keep_filtered)
    SNPEFF_ANNOTATE(PREPARE_VCF.out.vcf, SNPEFF_BUILD.out.db, params.genome_name, params.cohort_name)
    SNPSIFT_EXTRACT(SNPEFF_ANNOTATE.out.vcf, params.cohort_name)
    ANNOTATION_QC_REPORT(
        SNPEFF_ANNOTATE.out.vcf,
        PREPARE_VCF.out.stats,
        params.cohort_name,
        params.genome_name,
        vcf_path
    )
    if (params.ems) {
        EMS_REPORT(
            SNPEFF_ANNOTATE.out.vcf,
            optional_file(params.ems_alignment_qc),
            optional_file(params.ems_candidate_genes),
            optional_file(params.ems_homeolog_pairs),
            params.cohort_name
        )
    }
}