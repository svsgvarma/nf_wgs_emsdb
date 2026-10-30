process SNPEFF_BUILD {
    tag "${genome}"
    storeDir params.snpeff_db_store
    cpus 2
    memory '16 GB'
    input:
        path fasta
        path gff
        val genome
    output:
        path "${genome}", emit: db
    script:
    def xmx = Math.max(2, (task.memory.toGiga() as int) - 2)
    def is_gtf = gff.name ==~ /.*\.gtf(\.gz)?$/
    def annot_file = is_gtf ? 'genes.gtf' : 'genes.gff'
    def format = is_gtf ? '-gtf22' : '-gff3'
    """
    set -o pipefail
    SNPEFF_JAR=\${SNPEFF_JAR:-\$(find /opt/conda/share /usr/local/share -name snpEff.jar 2>/dev/null | head -n 1)}
    SNPEFF_HOME=\${SNPEFF_HOME:-\$(dirname "\$SNPEFF_JAR")}
    export SNPEFF_JAR SNPEFF_HOME
    mkdir -p data/${genome}
    case "${fasta}" in *.gz) zcat ${fasta} > data/${genome}/sequences.fa ;; *) ln -s "\$(readlink -f ${fasta})" data/${genome}/sequences.fa ;; esac
    case "${gff}" in
        *.gtf|*.gtf.gz) case "${gff}" in *.gz) zcat ${gff} | tr -d '\\r' > data/${genome}/${annot_file} ;; *) tr -d '\\r' < ${gff} > data/${genome}/${annot_file} ;; esac ;;
        *) clean_gff.sh ${gff} > data/${genome}/${annot_file} ;;
    esac
    make_snpeff_config.sh ${genome} '${params.snpeff_codon_tables}' > snpEff.config
    java -Xmx${xmx}g -jar \$SNPEFF_JAR build ${format} -noCheckCds -noCheckProtein -v \\
        -c snpEff.config -dataDir \$PWD/data ${genome} > build.log 2>&1 || { tail -n 50 build.log; exit 1; }
    mkdir -p ${genome}
    mv data/${genome}/snpEffectPredictor.bin ${genome}/
    mv build.log ${genome}/${genome}.build.log
    """
}