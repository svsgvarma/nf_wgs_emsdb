process PREPARE_REFERENCE {
    tag "${fasta.name}"
    storeDir "${params.index_store}/${fasta.baseName}"
    cpus 2
    memory { 16.GB * task.attempt }
    input:
        path fasta
    output:
        path "${fasta.name}.fai", emit: fai
        path "${fasta.baseName}.dict", emit: dict
        path "${fasta.name}.{amb,ann,bwt,pac,sa}", emit: bwa
    script:
    """
    samtools faidx ${fasta}
    gatk CreateSequenceDictionary -R ${fasta} -O ${fasta.baseName}.dict
    bwa index -a bwtsw ${fasta}
    """
}

process SCATTER_INTERVALS {
    tag "${fai.name}"
    publishDir "${params.outdir}/intervals", mode: 'copy', overwrite: true, pattern: 'shards.tsv'
    cpus 1
    memory '1 GB'
    input:
        path fai
        path intervals, stageAs: 'user_intervals/*'
    output:
        path "shard_*.bed", emit: shards
        path "shards.tsv", emit: summary
    script:
    def user = intervals.toString().endsWith('NO_FILE') ? '' : "--intervals '${intervals}'"
    def regex = params.contig_regex ? "--contig-regex '${params.contig_regex}'" : ''
    """
    make_shards.py --fai ${fai} ${user} --scatter-count ${params.scatter_count} ${regex} \\
        --min-contig-length ${params.min_contig_length} --min-fraction ${params.min_interval_fraction} --outdir .
    """
}