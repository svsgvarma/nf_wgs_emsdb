#!/usr/bin/env python3
"""Build a portable HTML and TSV summary for SNP-calling alignment and variant QC."""

import argparse
import glob
import gzip
import html
import os
from collections import Counter


def open_text(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, encoding='utf-8')


def flagstat_metrics(path):
    values = {}
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            columns = line.rstrip().split('\t')
            if len(columns) >= 3:
                values[columns[2]] = int(columns[0])
    primary = values.get('primary', 0)
    mapped = values.get('primary mapped', 0)
    duplicates = values.get('primary duplicates', values.get('duplicates', 0))
    return primary, mapped, duplicates


def samples(qc_dir):
    suffix = '.qc.flagstat.tsv'
    return [(os.path.basename(path)[:-len(suffix)], path)
            for path in sorted(glob.glob(os.path.join(qc_dir, '*' + suffix)))]


def variant_counts(path):
    counts = Counter()
    if not path:
        return counts
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            columns = line.rstrip().split('\t')
            if len(columns) < 5:
                continue
            ref, alts = columns[3], columns[4].split(',')
            counts['sites'] += 1
            if all(len(ref) == 1 and len(alt) == 1 for alt in alts):
                counts['snps'] += 1
            else:
                counts['indels'] += 1
    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--qc-dir', default='.')
    parser.add_argument('--vcf')
    parser.add_argument('--cohort', default='cohort')
    parser.add_argument('--outdir', default='.')
    parser.add_argument('--reference', default='')
    parser.add_argument('--intervals', default='')
    parser.add_argument('--min-mapping-rate', type=float, default=90)
    parser.add_argument('--min-mean-depth', type=float, default=10)
    parser.add_argument('--max-duplicate-rate', type=float, default=30)
    parser.add_argument('--max-missing-rate', type=float, default=10)
    parser.add_argument('--max-contigs', type=int, default=40)
    parser.add_argument('--chrom-pattern', default='')
    parser.add_argument('--many-samples', type=int, default=24)
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rows = []
    for sample, path in samples(args.qc_dir):
        primary, mapped, duplicates = flagstat_metrics(path)
        mapping_rate = 100 * mapped / primary if primary else 0
        duplicate_rate = 100 * duplicates / mapped if mapped else 0
        status = 'PASS' if mapping_rate >= args.min_mapping_rate and duplicate_rate <= args.max_duplicate_rate else 'FLAG'
        rows.append((sample, primary, mapped, mapping_rate, duplicates, duplicate_rate, status))
    prefix = os.path.join(args.outdir, args.cohort)
    with open(prefix + '.sample_alignment_qc.tsv', 'w', encoding='utf-8') as handle:
        handle.write('sample\tprimary_reads\tmapped_reads\tmapping_rate\tduplicate_reads\tduplicate_rate\tstatus\n')
        for row in rows:
            handle.write('\t'.join(str(value) for value in row) + '\n')
    variants = variant_counts(args.vcf)
    with open(prefix + '.variant_summary.tsv', 'w', encoding='utf-8') as handle:
        handle.write('metric\tvalue\n')
        for key in ('sites', 'snps', 'indels'):
            handle.write(f'{key}\t{variants[key]}\n')
    table_rows = ''.join('<tr>' + ''.join(f'<td>{html.escape(str(value))}</td>' for value in row) + '</tr>'
                         for row in rows)
    document = f'''<!doctype html><html><head><meta charset="utf-8"><title>SNP calling QC</title>
<style>body{{font-family:sans-serif;max-width:1100px;margin:2rem auto}}table{{border-collapse:collapse}}
th,td{{border:1px solid #ccc;padding:.4rem}}</style></head><body><h1>SNP-calling QC: {html.escape(args.cohort)}</h1>
<p>Reference: {html.escape(args.reference)}; intervals: {html.escape(args.intervals)}</p>
<p>Variant sites: {variants['sites']:,}; SNPs: {variants['snps']:,}; indels: {variants['indels']:,}</p>
<table><thead><tr><th>Sample</th><th>Primary</th><th>Mapped</th><th>Mapping %</th><th>Duplicates</th>
<th>Duplicate %</th><th>Status</th></tr></thead><tbody>{table_rows}</tbody></table></body></html>'''
    with open(prefix + '.snp_calling_qc_report.html', 'w', encoding='utf-8') as handle:
        handle.write(document)


if __name__ == '__main__':
    main()