/*
 * Call variants with GATK HaplotypeCaller
 */
process GATK_HAPLOTYPECALLER {

    //container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    // mode: 'symlink' does not work with S3 — symlinks are local filesystem constructs and cannot be created in S3. Using copy instead
    // publishDir params.outdir, mode: 'symlink'
    publishDir params.outdir, mode: 'copy'

    memory '16 GB'

    input:
        tuple path(input_bam), path(input_bam_index)
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list

    output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx

    script:
    """
    gatk HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -O ${input_bam}.g.vcf \
        -L ${interval_list} \
        -ERC GVCF
    """
}
