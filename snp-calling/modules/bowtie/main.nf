/*
 * Bowtie2 alignment process
 */

process BOWTIE_BUILD_INDEX {
    tag "$ref_fasta"
    publishDir params.outdir, mode: 'copy'

    input:
    path ref_fasta

    output:
    tuple val("bowtie_index"), path("bt2_index.*")

    script:
    """
    # Only build index if it does not exist
    if ls bt2_index.* 1> /dev/null 2>&1; then
        echo "Bowtie2 index already exists, skipping build."
    else
        bowtie2-build $ref_fasta bt2_index
    fi
    """
}

// Bowtie2 alignment process
process BOWTIE_ALIGN {
    tag "$sample_id"
    publishDir params.outdir, mode: 'symlink'

    input:
    tuple val(sample_id), path(reads1), path(reads2)
    tuple val(index_tag), path(index_files)

    output:
    tuple val(sample_id), path("${sample_id}.sam")

    script:
    """
    bowtie2 -x bt2_index -1 $reads1 -2 $reads2 -S ${sample_id}.sam
    """
}
