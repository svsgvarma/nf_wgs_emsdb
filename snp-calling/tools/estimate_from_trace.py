#!/usr/bin/env python3
import argparse
import csv
import re
import sys
from collections import defaultdict


def duration_hours(text):
    total = 0.0
    for value, unit in re.findall(r'([\d.]+)\s*(ms|d|h|m|s)', text or ''):
        value = float(value)
        total += {'ms': value / 3.6e6, 's': value / 3600, 'm': value / 60, 'h': value, 'd': value * 24}[unit]
    return total


def memory_gb(text):
    match = re.match(r'([\d.]+)\s*([KMGT]?B)', (text or '').strip())
    if not match:
        return 0.0
    return float(match.group(1)) * {'B': 1 / 1024 ** 3, 'KB': 1 / 1024 ** 2, 'MB': 1 / 1024,
                                    'GB': 1, 'TB': 1024}[match.group(2)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('trace')
    parser.add_argument('--target-samples', type=int, default=192)
    parser.add_argument('--size-factor', type=float, default=1.0)
    parser.add_argument('--gb-per-vcpu', type=float, default=4.0)
    parser.add_argument('--spot-price', type=float, default=0.018)
    parser.add_argument('--ondemand-price', type=float, default=0.048)
    args = parser.parse_args()
    with open(args.trace, encoding='utf-8') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if not rows:
        sys.exit('empty trace')
    completed = [row for row in rows if row.get('status') in ('COMPLETED', 'CACHED')]
    pilot_samples = sum(1 for row in completed if 'BWA_ALIGN' in (row.get('process') or row.get('name', ''))) or 1
    totals = defaultdict(float)
    for row in completed:
        process = (row.get('process') or row.get('name', '')).split(':')[-1].split(' (')[0]
        billed = max(float(row.get('cpus') or 1), memory_gb(row.get('memory')) / args.gb_per_vcpu)
        totals[process] += billed * duration_hours(row.get('realtime'))
    factor = args.target_samples / pilot_samples * args.size_factor
    target = sum(totals.values()) * factor
    print(f'Pilot samples: {pilot_samples}')
    print(f'Target billed compute: {target:,.1f} vCPU-hours')
    print(f'Spot estimate: ${target * args.spot_price:,.0f}')
    print(f'On-demand estimate: ${target * args.ondemand_price:,.0f}')


if __name__ == '__main__':
    main()