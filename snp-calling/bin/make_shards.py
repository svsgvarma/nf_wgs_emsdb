#!/usr/bin/env python3
import argparse
import os
import re
import sys


def read_fai(path):
    contigs = []
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            columns = line.rstrip('\r\n').split('\t')
            if len(columns) >= 2 and columns[0]:
                contigs.append((columns[0], int(columns[1])))
    if not contigs:
        sys.exit(f'ERROR: no contigs found in {path}')
    return contigs


def read_intervals(path):
    name = path.lower()
    intervals = []
    with open(path, encoding='utf-8') as handle:
        for number, line in enumerate(handle, 1):
            line = line.rstrip('\r\n')
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            if name.endswith('.interval_list'):
                if line.startswith('@'):
                    continue
                columns = line.split('\t')
                intervals.append((columns[0], int(columns[1]) - 1, int(columns[2])))
            elif name.endswith(('.list', '.intervals')):
                match = re.match(r'^([^:\s]+)(?::(\d[\d,]*)(?:-(\d[\d,]*))?)?$', line.strip())
                if not match:
                    sys.exit(f'ERROR: {os.path.basename(path)} line {number}: cannot parse interval "{line}"')
                contig, start, end = match.group(1), match.group(2), match.group(3)
                if start is None:
                    intervals.append((contig, 0, None))
                else:
                    start = int(start.replace(',', ''))
                    end = int(end.replace(',', '')) if end else start
                    intervals.append((contig, start - 1, end))
            else:
                columns = line.split('\t') if '\t' in line else line.split()
                if len(columns) < 3:
                    sys.exit(f'ERROR: {os.path.basename(path)} line {number}: a BED line needs contig, start and end')
                intervals.append((columns[0], int(columns[1]), int(columns[2])))
    return intervals


def merge(intervals):
    merged = []
    for contig, start, end in intervals:
        if merged and merged[-1][0] == contig and start <= merged[-1][2]:
            merged[-1] = (contig, merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append((contig, start, end))
    return merged


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fai', required=True)
    parser.add_argument('--intervals', default='')
    parser.add_argument('--scatter-count', type=int, default=12)
    parser.add_argument('--contig-regex', default='')
    parser.add_argument('--min-contig-length', type=int, default=0)
    parser.add_argument('--min-fraction', type=float, default=0.01)
    parser.add_argument('--outdir', default='.')
    args = parser.parse_args()
    contigs = read_fai(args.fai)
    order = {contig: index for index, (contig, _) in enumerate(contigs)}
    lengths = dict(contigs)
    genome_size = sum(lengths.values())
    if args.intervals:
        raw = read_intervals(args.intervals)
        missing = sorted({contig for contig, _, _ in raw if contig not in lengths})
        if missing:
            sys.exit(f'ERROR: intervals name contigs absent from the reference: {", ".join(missing[:10])}')
        regions = [(contig, max(0, start), lengths[contig] if end is None else min(end, lengths[contig]))
                   for contig, start, end in raw]
        regions = merge(sorted((region for region in regions if region[2] > region[1]),
                               key=lambda region: (order[region[0]], region[1], region[2])))
    else:
        pattern = re.compile(args.contig_regex) if args.contig_regex else None
        regions = [(contig, 0, length) for contig, length in contigs
                   if length >= args.min_contig_length and (pattern is None or pattern.search(contig))]
    total = sum(end - start for _, start, end in regions)
    if not regions or not total:
        sys.exit('ERROR: no calling regions left')
    fraction = total / genome_size
    if args.min_fraction and fraction < args.min_fraction:
        sys.exit(f'ERROR: calling regions cover only {100 * fraction:.3f}% of the reference')
    target = total / max(1, args.scatter_count)
    shards, current, current_size = [], [], 0
    for region in regions:
        current.append(region)
        current_size += region[2] - region[1]
        if current_size >= target and len(shards) < max(1, args.scatter_count) - 1:
            shards.append(current)
            current, current_size = [], 0
    if current:
        shards.append(current)
    os.makedirs(args.outdir, exist_ok=True)
    width = max(3, len(str(len(shards) - 1)))
    with open(os.path.join(args.outdir, 'shards.tsv'), 'w', encoding='utf-8') as summary:
        summary.write('shard\tintervals\tcontigs\tbp\twhole_contigs\n')
        for index, shard in enumerate(shards):
            name = f'shard_{index:0{width}d}'
            with open(os.path.join(args.outdir, f'{name}.bed'), 'w', encoding='utf-8') as handle:
                for contig, start, end in shard:
                    handle.write(f'{contig}\t{start}\t{end}\n')
            whole = all(start == 0 and end == lengths[contig] for contig, start, end in shard)
            summary.write(f'{name}\t{len(shard)}\t{len({r[0] for r in shard})}\t'
                          f'{sum(r[2] - r[1] for r in shard)}\t{str(whole).lower()}\n')


if __name__ == '__main__':
    main()