process SNPEFF_ANNOTATE {
    tag "${cohort}"
    publishDir "${params.outdir}/annotation", mode: 'copy', overwrite: true
    cpus 4
    memory '16 GB'
    input:
        tuple path(vcf), path(vcf_index)
        path db
        val genome
        val cohort
    output:
        tuple path("${cohort}.snpeff.vcf.gz"), path("${cohort}.snpeff.vcf.gz.tbi"), emit: vcf
        path "${cohort}.snpeff_summary.*", emit: summary
    script:
    def xmx = Math.max(2, (task.memory.toGiga() as int) - 2)
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
        ${params.snpeff_args} ${genome} ${vcf} | bgzip -@ ${task.cpus} > ${cohort}.snpeff.vcf.gz
    tabix -p vcf ${cohort}.snpeff.vcf.gz
    """
}