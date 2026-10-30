#!/usr/bin/env python3
"""Generate general annotation QC outputs from a SnpEff-annotated VCF."""

import argparse
import gzip
import html
import os
from collections import Counter

from report_lib import ann_record, first_ann, sub_class, write_tsv


def parse_vcf(path):
    opener = gzip.open if path.endswith('.gz') else open
    counts, impacts, genes, chromosomes = Counter(), Counter(), Counter(), Counter()
    samples = []
    with opener(path, 'rt', encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('#CHROM'):
                samples = line.rstrip().split('\t')[9:]
                continue
            if line.startswith('#'):
                continue
            columns = line.rstrip().split('\t')
            if len(columns) < 8 or columns[6] not in ('PASS', '.'):
                continue
            counts['variants'] += 1
            change = sub_class(columns[3], columns[4].split(',')[0])
            if change:
                counts['snps'] += 1
            fields = first_ann(columns[7])
            if fields:
                annotation = ann_record(fields)
                impacts[annotation['impact']] += 1
                if annotation['gene']:
                    genes[(annotation['gene'], annotation['impact'])] += 1
            chromosomes[columns[0]] += 1
    return counts, impacts, genes, chromosomes, samples


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', required=True)
    parser.add_argument('--stats-dir', default='.')
    parser.add_argument('--cohort', default='cohort')
    parser.add_argument('--genome', default='')
    parser.add_argument('--input-vcf', default='')
    parser.add_argument('--snpeff-summary', default='')
    parser.add_argument('--max-contigs', type=int, default=40)
    parser.add_argument('--many-samples', type=int, default=24)
    parser.add_argument('--top-genes', type=int, default=50)
    parser.add_argument('--window', type=int, default=1000000)
    parser.add_argument('--outdir', default='.')
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    counts, impacts, genes, chromosomes, samples = parse_vcf(args.vcf)
    prefix = os.path.join(args.outdir, args.cohort)
    write_tsv(prefix + '.annotation_summary.tsv', ['metric', 'value'], list(counts.items()) + list(impacts.items()))
    write_tsv(prefix + '.impact_by_chromosome.tsv', ['chromosome', 'variants'], chromosomes.items())
    write_tsv(prefix + '.genes_high_moderate.tsv', ['gene', 'impact', 'variants'],
              [(gene, impact, count) for (gene, impact), count in genes.items() if impact in ('HIGH', 'MODERATE')])
    write_tsv(prefix + '.sample_summary.tsv', ['sample'], [(sample,) for sample in samples])
    write_tsv(prefix + '.variant_density.tsv', ['chromosome', 'window', 'variants'],
              [(chromosome, 0, count) for chromosome, count in chromosomes.items()])
    rows = ''.join(f'<tr><td>{html.escape(key)}</td><td>{value}</td></tr>' for key, value in counts.items())
    with open(prefix + '.annotation_qc_report.html', 'w', encoding='utf-8') as handle:
        handle.write(f'<!doctype html><html><body><h1>Annotation QC: {html.escape(args.cohort)}</h1>'
                     f'<p>Genome: {html.escape(args.genome)}</p><table>{rows}</table></body></html>')


if __name__ == '__main__':
    main()