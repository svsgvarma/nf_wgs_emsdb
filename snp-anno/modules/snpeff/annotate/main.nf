process SNPEFF_ANNOTATE {
    tag "${cohort}"
    publishDir "${params.outdir}/annotation", mode: 'copy', overwrite: true, pattern: '*.snpeff_summary.*'
    cpus 4
    memory '16 GB'
    input:
        tuple path(vcf), path(vcf_index)
        path db
        val genome
        val cohort
    output:
        path "${cohort}.snpeff.vcf", emit: vcf
        path "${cohort}.snpeff_summary.*", emit: summary
    script:
    def xmx = Math.max(2, (task.memory.toGiga() as int) - 2)
    // Plain VCF out: the SnpEff image has no bgzip/tabix, so BGZIP_INDEX_VCF compresses it.
    """
    set -o pipefail
    SNPEFF_JAR=\${SNPEFF_JAR:-\$(find /opt/conda/share /usr/local/share -name snpEff.jar 2>/dev/null | head -n 1)}
    SNPEFF_HOME=\${SNPEFF_HOME:-\$(dirname "\$SNPEFF_JAR")}
    export SNPEFF_JAR SNPEFF_HOME
    make_snpeff_config.sh ${genome} '${params.snpeff_codon_tables}' > snpEff.config
    mkdir -p data
    ln -s "\$(readlink -f ${db})" data/${genome}
    java -Xmx${xmx}g -jar \$SNPEFF_JAR ann -nodownload -c snpEff.config -dataDir \$PWD/data \\
        -stats ${cohort}.snpeff_summary.html -csvStats ${cohort}.snpeff_summary.csv \\
        ${params.snpeff_args} ${genome} ${vcf} > ${cohort}.snpeff.vcf
    """
}

process BGZIP_INDEX_VCF {
    tag "${vcf.name}"
    publishDir "${params.outdir}/annotation", mode: 'copy', overwrite: true
    cpus 2
    memory '2 GB'
    input:
        path vcf
    output:
        tuple path("${vcf.name}.gz"), path("${vcf.name}.gz.tbi"), emit: vcf
    script:
    """
    bgzip -@ ${task.cpus} -c ${vcf} > ${vcf.name}.gz
    tabix -p vcf ${vcf.name}.gz
    """
}
