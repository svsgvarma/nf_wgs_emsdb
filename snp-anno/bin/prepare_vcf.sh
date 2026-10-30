#!/usr/bin/env bash
set -euo pipefail
vcf="$1"; fasta="$2"; cohort="$3"; normalize="$4"; hard_filter="$5"; keep_filtered="$6"; threads="${7:-2}"
bcftools stats --threads "$threads" "$vcf" > "${cohort}.prepare_1_input.stats.txt"
if [[ "$normalize" == "true" ]]; then
    case "$fasta" in *.gz) zcat "$fasta" > ref.fa ;; *) ln -sf "$(readlink -f "$fasta")" ref.fa ;; esac
    samtools faidx ref.fa
    bcftools norm --threads "$threads" -m -any --check-ref w -f ref.fa -Oz -o norm.vcf.gz "$vcf"
else
    bcftools view --threads "$threads" -Oz -o norm.vcf.gz "$vcf"
fi
tabix -f -p vcf norm.vcf.gz
bcftools stats --threads "$threads" norm.vcf.gz > "${cohort}.prepare_2_normalized.stats.txt"
bcftools view -h norm.vcf.gz > header.txt
build_expr() {
    local expr="" cond tag
    for cond in "$@"; do
        tag="${cond%% *}"
        if [[ "$tag" == "QUAL" ]]; then expr="${expr:+$expr || }$cond"
        elif grep -q "^##INFO=<ID=${tag}," header.txt; then expr="${expr:+$expr || }INFO/$cond"; fi
    done
    echo "$expr"
}
if [[ "$hard_filter" == "true" ]]; then
    snp_expr=$(build_expr "QD < 2.0" "QUAL < 30" "SOR > 3.0" "FS > 60.0" "MQ < 40.0" "MQRankSum < -12.5" "ReadPosRankSum < -8.0")
    indel_expr=$(build_expr "QD < 2.0" "QUAL < 30" "FS > 200.0" "SOR > 10.0" "ReadPosRankSum < -20.0")
    bcftools filter --threads "$threads" -m+ -s SnpHardFilter -e "TYPE=\"snp\" && ($snp_expr)" -Ou norm.vcf.gz \
        | bcftools filter --threads "$threads" -m+ -s IndelHardFilter -e "TYPE=\"indel\" && ($indel_expr)" -Oz -o soft.vcf.gz
    bcftools query -f '%FILTER\n' soft.vcf.gz | awk '{c[$1]++} END {for (k in c) print k "\t" c[k]}' > "${cohort}.prepare_filter_counts.tsv"
    if [[ "$keep_filtered" == "true" ]]; then mv soft.vcf.gz "${cohort}.prepared.vcf.gz"
    else bcftools view --threads "$threads" -f PASS,. -Oz -o "${cohort}.prepared.vcf.gz" soft.vcf.gz; rm -f soft.vcf.gz; fi
else
    : > "${cohort}.prepare_filter_counts.tsv"
    mv norm.vcf.gz "${cohort}.prepared.vcf.gz"
fi
rm -f norm.vcf.gz norm.vcf.gz.tbi header.txt
tabix -f -p vcf "${cohort}.prepared.vcf.gz"
bcftools stats --threads "$threads" "${cohort}.prepared.vcf.gz" > "${cohort}.prepare_3_prepared.stats.txt"