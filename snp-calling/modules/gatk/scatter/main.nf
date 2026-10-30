def java_heap(task, reserve_gb) {
    return Math.max(1, (task.memory.toGiga() as int) - reserve_gb)
}

process GATK_HAPLOTYPECALLER_SHARD {
    tag "${sample_id}:${shard_id}"
    cpus 2
    memory { 6.GB * task.attempt }
    input:
        tuple val(sample_id), path(bam), path(bai), val(shard_id), path(shard_bed)
        path ref_fasta
        path ref_index
        path ref_dict
    output:
        tuple val(shard_id), val(sample_id), path("${sample_id}.${shard_id}.g.vcf.gz"), path("${sample_id}.${shard_id}.g.vcf.gz.tbi"), emit: gvcf
    script:
    """
    gatk --java-options "-Xmx${java_heap(task, 2)}g" HaplotypeCaller \\
        -R ${ref_fasta} -I ${bam} -L ${shard_bed} -O ${sample_id}.${shard_id}.g.vcf.gz \\
        -ERC GVCF -ploidy ${params.ploidy} --native-pair-hmm-threads ${task.cpus} \\
        --tmp-dir . ${params.haplotypecaller_args}
    """
}

process GATK_JOINT_GENOTYPE_SHARD {
    tag "${shard_id}"
    cpus 4
    memory { 16.GB * task.attempt }
    input:
        tuple val(shard_id), val(sample_ids), path(gvcfs), path(tbis), path(shard_bed)
        path ref_fasta
        path ref_index
        path ref_dict
        val cohort_name
    output:
        tuple val(shard_id), path("${cohort_name}.${shard_id}.vcf.gz"), path("${cohort_name}.${shard_id}.vcf.gz.tbi"), emit: vcf
    script:
    def heap = Math.max(2, (task.memory.toGiga() * 0.6) as int)
    def gvcf_list = gvcfs instanceof List ? gvcfs : [gvcfs]
    def id_list = sample_ids instanceof List ? sample_ids : [sample_ids]
    def sample_map = [id_list, gvcf_list].transpose().collect { sample, gvcf -> "${sample}\t${gvcf}" }.join('\n')
    """
    printf '%s\\n' "${sample_map}" > sample_map.tsv
    mkdir -p tmp
    MERGE=""
    N_INTERVALS=\$(wc -l < ${shard_bed})
    WHOLE=\$(awk 'NR==FNR { len[\$1] = \$2; next } !(\$2 == 0 && \$3 == len[\$1]) { bad = 1 } END { print bad ? "no" : "yes" }' ${ref_index} ${shard_bed})
    if [ "\$WHOLE" = "yes" ] && [ "\$N_INTERVALS" -gt ${params.genomicsdb_merge_contigs_above} ]; then
        MERGE="--merge-contigs-into-num-partitions 25"
    fi
    gatk --java-options "-Xmx${heap}g" GenomicsDBImport \\
        --sample-name-map sample_map.tsv --genomicsdb-workspace-path ${shard_id}_gdb \\
        -L ${shard_bed} --batch-size ${params.genomicsdb_batch_size} \\
        --reader-threads ${task.cpus} --tmp-dir tmp \$MERGE
    gatk --java-options "-Xmx${heap}g" GenotypeGVCFs \\
        -R ${ref_fasta} -V gendb://${shard_id}_gdb -L ${shard_bed} \\
        -O ${cohort_name}.${shard_id}.vcf.gz -stand-call-conf ${params.genotype_call_conf} --tmp-dir tmp
    rm -rf ${shard_id}_gdb tmp
    """
}

process GATK_GATHER_VCFS {
    tag "${out_name}"
    publishDir { "${params.outdir}${publish_subdir ? '/' + publish_subdir : ''}" }, mode: 'copy', overwrite: true
    cpus 1
    memory { 4.GB * task.attempt }
    input:
        tuple val(out_name), path(vcfs), path(tbis)
        val publish_subdir
    output:
        tuple val(out_name), path("${out_name}"), path("${out_name}.tbi"), emit: vcf
    script:
    def inputs = (vcfs instanceof List ? vcfs : [vcfs]).collect { "-I ${it}" }.join(' ')
    """
    gatk --java-options "-Xmx${java_heap(task, 1)}g" GatherVcfs ${inputs} -O ${out_name}
    gatk IndexFeatureFile -I ${out_name}
    """
}