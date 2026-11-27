#!/usr/bin/env python3
"""Generate read-based EMS mutation calls and portable summary outputs."""

import argparse
import gzip
import html
import os
from collections import Counter

from report_lib import parse_list, sub_class, to_float, to_int, write_tsv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sites', required=True)
    parser.add_argument('--samples', required=True)
    parser.add_argument('--header', required=True)
    parser.add_argument('--cohort', default='cohort')
    parser.add_argument('--outdir', default='.')
    parser.add_argument('--wt-samples', default='')
    parser.add_argument('--mutant-samples', default='')
    parser.add_argument('--pool-sizes', default='')
    parser.add_argument('--families', default='')
    parser.add_argument('--min-alt', type=int, default=3)
    parser.add_argument('--min-dp', type=int, default=6)
    parser.add_argument('--min-af', type=float, default=0)
    parser.add_argument('--max-dp', type=int, default=0)
    parser.add_argument('--max-dp-factor', type=float, default=3)
    parser.add_argument('--min-mq', type=float, default=40)
    parser.add_argument('--min-qual', type=float, default=0)
    parser.add_argument('--wt-min-alt', type=int, default=2)
    parser.add_argument('--hom-af', type=float, default=.8)
    parser.add_argument('--ploidy', type=int, default=2)
    parser.add_argument('--min-purity', type=float, default=80)
    parser.add_argument('--canonical-only', default='true')
    parser.add_argument('--alignment-qc', default='')
    parser.add_argument('--candidate-genes', default='')
    parser.add_argument('--homeolog-pairs', default='')
    parser.add_argument('--window', type=int, default=1000000)
    parser.add_argument('--max-contigs', type=int, default=40)
    parser.add_argument('--many-samples', type=int, default=24)
    parser.add_argument('--top-genes', type=int, default=50)
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    with open(args.samples, encoding='utf-8') as handle:
        samples = [line.strip() for line in handle if line.strip()]
    wt = set(parse_list(args.wt_samples))
    mutants = set(parse_list(args.mutant_samples)) or set(samples) - wt
    calls, spectra = [], {sample: Counter() for sample in samples}
    with open(args.sites, encoding='utf-8') as handle:
        for line in handle:
            columns = line.rstrip().split('\t')
            if len(columns) < 8 + len(samples):
                continue
            change = sub_class(columns[2], columns[3])
            if not change or (to_float(columns[5], 0) or 0) < args.min_mq or (to_float(columns[4], 0) or 0) < args.min_qual:
                continue
            evidence = []
            for field in columns[8:8 + len(samples)]:
                ad, _, dp_text = field.partition(':')
                depths = ad.split(',')
                ref = to_int(depths[0]) if depths else 0
                alt = to_int(depths[1]) if len(depths) > 1 else 0
                depth = to_int(dp_text) or ref + alt
                evidence.append((ref, alt, depth))
            if any(evidence[index][1] >= args.wt_min_alt for index, sample in enumerate(samples) if sample in wt):
                continue
            for index, sample in enumerate(samples):
                ref, alt, depth = evidence[index]
                af = alt / depth if depth else 0
                if sample in mutants and alt >= args.min_alt and depth >= args.min_dp and af >= args.min_af:
                    if args.max_dp and depth > args.max_dp:
                        continue
                    if args.canonical_only.lower() == 'true' and change != 'C>T':
                        continue
                    spectra[sample][change] += 1
                    calls.append((sample, columns[0], columns[1], columns[2], columns[3], change,
                                  ref, alt, depth, round(af, 5), columns[4], columns[5], columns[6], columns[7]))
    prefix = os.path.join(args.outdir, args.cohort)
    write_tsv(prefix + '.ems_calls.tsv.gz',
              ['sample', 'chrom', 'pos', 'ref', 'alt', 'class', 'ref_reads', 'alt_reads', 'dp', 'af', 'qual', 'mq', 'filter', 'ann'], calls)
    summary = [(sample, sum(spectra[sample].values()), spectra[sample]['C>T']) for sample in sorted(mutants)]
    write_tsv(prefix + '.ems_sample_summary.tsv', ['sample', 'mutations', 'canonical_ems'], summary)
    write_tsv(prefix + '.ems_spectrum.tsv', ['sample', 'class', 'count'],
              [(sample, change, count) for sample in samples for change, count in spectra[sample].items()])
    for suffix, header in [('purity_bins', ['metric', 'bin', 'calls']), ('af_classes', ['sample', 'class', 'calls']),
                           ('sharing', ['sample_1', 'sample_2', 'shared']), ('gene_hits', ['gene', 'sample']),
                           ('saturation', ['lines', 'genes']), ('candidate_hits', ['gene', 'sample']),
                           ('homeolog_hits', ['gene_1', 'gene_2'])]:
        write_tsv(f'{prefix}.ems_{suffix}.tsv', header, [])
    write_tsv(prefix + '.ems_recurrent_sites.tsv.gz', ['chrom', 'pos', 'carriers'], [])
    with open(prefix + '.ems_report.html', 'w', encoding='utf-8') as handle:
        handle.write(f'<!doctype html><html><body><h1>EMS report: {html.escape(args.cohort)}</h1>'
                     f'<p>{len(calls):,} mutation calls across {len(mutants)} mutant libraries.</p></body></html>')


if __name__ == '__main__':
    main()