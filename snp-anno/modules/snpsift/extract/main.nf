process SNPSIFT_EXTRACT {
    tag "${cohort}"
    publishDir "${params.outdir}/annotation", mode: 'copy', overwrite: true
    cpus 1
    memory '8 GB'
    input:
        tuple path(vcf), path(vcf_index)
        val cohort
    output:
        path "${cohort}.snpeff_effects.tsv.gz", emit: tsv
    script:
    def xmx = Math.max(2, (task.memory.toGiga() as int) - 2)
    def fields = params.snpsift_fields.trim().split(/\s+/).collect { "'${it}'" }.join(' ')
    """
    set -o pipefail
    SNPSIFT_JAR=\${SNPSIFT_JAR:-\$(find /opt/conda/share /usr/local/share -name SnpSift.jar 2>/dev/null | head -n 1)}
    export SNPSIFT_JAR
    java -Xmx${xmx}g -jar \$SNPSIFT_JAR extractFields -s "," -e "." ${vcf} ${fields} | gzip > ${cohort}.snpeff_effects.tsv.gz
    """
}