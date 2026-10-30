#!/usr/bin/env bash
set -euo pipefail
genome="$1"
codon_tables="${2:-}"
cat "${SNPEFF_HOME}/snpEff.config"
echo
echo "# Custom genome added by the snp-anno pipeline"
echo "${genome}.genome : ${genome}"
if [[ -n "${codon_tables}" ]]; then
    IFS=',' read -ra pairs <<< "${codon_tables}"
    for pair in "${pairs[@]}"; do
        chrom="${pair%%:*}"
        table="${pair#*:}"
        echo "${genome}.${chrom// /}.codonTable : ${table// /}"
    done
fi