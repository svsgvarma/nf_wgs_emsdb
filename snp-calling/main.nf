#!/usr/bin/env nextflow

params.reads_fastq = "${params.input_dir}/data/fastq/NIST_Illumina_2x250bps_*/*_R{1,2}_*.fastq.gz"
params.reference = "${params.input_dir}/data/ref/GRCh38_GIAB/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
params.reference_index = null
params.reference_index_bwa = null
params.reference_dict = null
params.intervals = null
params.contig_regex = ''
params.min_contig_length = 0
params.min_interval_fraction = 0.01
params.scatter_count = 12
params.ploidy = 2
params.haplotypecaller_args = ''
params.genotype_call_conf = 30
params.genomicsdb_batch_size = 50
params.genomicsdb_merge_contigs_above = 100
params.gather_gvcfs = true
params.cohort_name = "family_trio"
params.publish_bams = false
params.qc_min_mapping_rate = 90
params.qc_min_mean_depth = 10
params.qc_max_duplicate_rate = 30
params.qc_max_missing_rate = 10
params.qc_max_contigs = 40
params.qc_chrom_pattern = ''
params.qc_many_samples = 24
params.qc_report_cpus = 8

include { BWA_INDEX; BWA_ALIGN } from './modules/bwa/main.nf'
include { SAMTOOLS_SORT; SAMTOOLS_ADDRG; SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'
include { SAMTOOLS_MARKDUP } from './modules/samtools/markdup/main.nf'
include { SAMTOOLS_QC } from './modules/samtools/qc/main.nf'
include { PREPARE_REFERENCE; SCATTER_INTERVALS } from './modules/reference/main.nf'
include { GATK_HAPLOTYPECALLER_SHARD; GATK_JOINT_GENOTYPE_SHARD } from './modules/gatk/scatter/main.nf'
include { GATK_GATHER_VCFS as GATHER_SAMPLE_GVCF; GATK_GATHER_VCFS as GATHER_JOINT_VCF } from './modules/gatk/scatter/main.nf'
include { SNP_QC_REPORT } from './modules/qc_report/main.nf'
include { MULTIQC } from './modules/multiqc/main.nf'

workflow {
    reads_ch = Channel.fromFilePairs(params.reads_fastq, flat: true)
        .ifEmpty { error "No paired FASTQ files match --reads_fastq ${params.reads_fastq}" }

    ref_file = file(params.reference, checkIfExists: true)
    if (ref_file.name.endsWith('.gz'))
        error "--reference must be an uncompressed FASTA (GATK needs random access): ${params.reference}"
    def ref_stem = params.reference.replaceFirst(/\.(fa|fasta|fna)$/, '')
    def fai_path = params.reference_index ?: "${params.reference}.fai"
    def dict_path = params.reference_dict ?: "${ref_stem}.dict"
    def bwa_prefix = params.reference_index_bwa ?: params.reference
    def bwa_ext = ['amb', 'ann', 'bwt', 'pac', 'sa']
    def have_all = file(fai_path).exists() && file(dict_path).exists() &&
                     bwa_ext.every { ext -> file("${bwa_prefix}.${ext}").exists() }

    if (have_all) {
        ref_index_ch = Channel.value(file(fai_path))
        ref_dict_ch = Channel.value(file(dict_path))
        bwa_index_ch = Channel.value(bwa_ext.collect { ext -> file("${bwa_prefix}.${ext}") })
    } else {
        log.info "Reference indexes not found next to ${ref_file.name}; building them in ${params.index_store}/${ref_file.baseName}"
        PREPARE_REFERENCE(ref_file)
        ref_index_ch = PREPARE_REFERENCE.out.fai
        ref_dict_ch = PREPARE_REFERENCE.out.dict
        bwa_index_ch = PREPARE_REFERENCE.out.bwa
    }

    intervals_in = params.intervals ? file(params.intervals, checkIfExists: true) : file("${projectDir}/assets/NO_FILE")
    SCATTER_INTERVALS(ref_index_ch, intervals_in)
    shards_ch = SCATTER_INTERVALS.out.shards.flatten().map { bed -> tuple(bed.baseName, bed) }
    n_shards_ch = shards_ch.count()
    n_samples_ch = reads_ch.count()

    bwa_out = BWA_ALIGN(reads_ch, bwa_index_ch)
    SAMTOOLS_MARKDUP(bwa_out)
    indexed_bam_ch = SAMTOOLS_INDEX(SAMTOOLS_MARKDUP.out.bam)
    SAMTOOLS_QC(indexed_bam_ch)

    GATK_HAPLOTYPECALLER_SHARD(indexed_bam_ch.combine(shards_ch), ref_file, ref_index_ch, ref_dict_ch)

    joint_in = GATK_HAPLOTYPECALLER_SHARD.out.gvcf
        .combine(n_samples_ch)
        .map { shard_id, sample_id, gvcf, tbi, n -> tuple(groupKey(shard_id, n as int), sample_id, gvcf, tbi) }
        .groupTuple()
        .map { key, samples, gvcfs, tbis -> tuple(key.toString(), samples, gvcfs, tbis) }
        .join(shards_ch)
    GATK_JOINT_GENOTYPE_SHARD(joint_in, ref_file, ref_index_ch, ref_dict_ch, params.cohort_name)

    joint_gather_in = GATK_JOINT_GENOTYPE_SHARD.out.vcf
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { rows -> tuple("${params.cohort_name}.joint.vcf.gz", rows.collect { it[1] }, rows.collect { it[2] }) }
    GATHER_JOINT_VCF(joint_gather_in, '')

    if (params.gather_gvcfs) {
        sample_gather_in = GATK_HAPLOTYPECALLER_SHARD.out.gvcf
            .combine(n_shards_ch)
            .map { shard_id, sample_id, gvcf, tbi, n -> tuple(groupKey(sample_id, n as int), shard_id, gvcf, tbi) }
            .groupTuple()
            .map { key, shard_ids, gvcfs, tbis ->
                def rows = [shard_ids, gvcfs, tbis].transpose().sort { it[0] }
                tuple("${key}.g.vcf.gz", rows.collect { it[1] }, rows.collect { it[2] })
            }
        GATHER_SAMPLE_GVCF(sample_gather_in, 'gvcf')
    }

    SNP_QC_REPORT(
        SAMTOOLS_QC.out.qc.mix(SAMTOOLS_MARKDUP.out.metrics).collect(),
        GATHER_JOINT_VCF.out.vcf.map { name, vcf, tbi -> vcf },
        params.cohort_name
    )
}