process SAMTOOLS_QC {
    tag "${sample_id}"
    cpus 4
    memory { 4.GB * task.attempt }
    input:
        tuple val(sample_id), path(bam_file), path(bam_index)
    output:
        path "${sample_id}.qc.*", emit: qc
    script:
    """
    samtools flagstat -@ ${task.cpus} -O tsv ${bam_file} > ${sample_id}.qc.flagstat.tsv
    samtools stats -@ ${task.cpus} -d ${bam_file} > ${sample_id}.qc.stats.txt
    samtools coverage ${bam_file} > ${sample_id}.qc.coverage.tsv
    samtools view -@ ${task.cpus} -c -f 0x2 -F 0xF0C -q 20 ${bam_file} > ${sample_id}.qc.hq_reads.txt
    """
}