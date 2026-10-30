"""Shared standard-library helpers for snp-anno HTML and TSV reports."""

import gzip
import html
import json
import math
from datetime import datetime

IMPACTS = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
IMPACT_COLORS = {'HIGH': '#c0392b', 'MODERATE': '#e67e22', 'LOW': '#f1c40f', 'MODIFIER': '#95a5a6'}
COLORS = ['#3498db', '#e67e22', '#2ecc71', '#9b59b6', '#e74c3c', '#1abc9c']
SUB_CLASSES = ['C>T', 'C>A', 'C>G', 'T>C', 'T>A', 'T>G']
SUB_COLORS = dict(zip(SUB_CLASSES, COLORS))
REGIONS = ['coding/splice', 'UTR', 'intron', 'flanking', 'intergenic', 'other']
REGION_COLORS = dict(zip(REGIONS, COLORS))
_COMP = str.maketrans('ACGT', 'TGCA')


def to_int(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def to_float(value, default=None):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def pct(value, total):
    return 100 * value / total if total else 0


def fmt_int(value):
    return f'{int(round(value)):,}'


def fmt_f(value, digits=2):
    return f'{value:,.{digits}f}'


def parse_list(text):
    return [value.strip() for value in (text or '').replace(' ', ',').split(',') if value.strip()]


def parse_pairs(text, cast=str):
    return {key.strip(): cast(value.strip()) for item in (text or '').replace(' ', ',').split(',') if ':' in item
            for key, value in [item.split(':', 1)]}


def write_tsv(path, header, rows):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt', encoding='utf-8') as handle:
        handle.write('\t'.join(header) + '\n')
        for row in rows:
            handle.write('\t'.join('' if value is None else str(value) for value in row) + '\n')


def sub_class(ref, alt):
    if len(ref) != 1 or len(alt) != 1 or ref not in 'ACGT' or alt not in 'ACGT' or ref == alt:
        return None
    if ref in 'AG':
        ref, alt = ref.translate(_COMP), alt.translate(_COMP)
    return f'{ref}>{alt}'


def is_transition(change):
    return change in ('C>T', 'T>C')


def first_ann(info):
    value = info.split('ANN=', 1)[1].split(';', 1)[0] if 'ANN=' in info else info
    fields = value.split(',', 1)[0].split('|')
    return fields if len(fields) >= 4 else None


def ann_record(fields):
    fields = fields + [''] * (16 - len(fields))
    return {'effect': fields[1], 'impact': fields[2], 'gene': fields[3] or fields[4], 'gene_id': fields[4],
            'feature': fields[6], 'biotype': fields[7], 'hgvs_c': fields[9], 'hgvs_p': fields[10],
            'messages': fields[15]}


def contig_lengths(lines):
    lengths = {}
    for line in lines:
        if line.startswith('##contig=<'):
            fields = dict(item.split('=', 1) for item in line.strip()[10:-1].split(',') if '=' in item)
            lengths[fields.get('ID', '')] = to_int(fields.get('length'))
    return lengths


def ad_dp(field, ad_index, dp_index):
    parts = field.split(':')
    depths = parts[ad_index].split(',') if 0 <= ad_index < len(parts) else []
    ref = to_int(depths[0]) if depths else 0
    alt = to_int(depths[1]) if len(depths) > 1 else 0
    depth = to_int(parts[dp_index]) if 0 <= dp_index < len(parts) else ref + alt
    return ref, alt, depth or ref + alt


def html_table(header, rows, table_id=None):
    headings = ''.join(f'<th>{html.escape(str(value))}</th>' for value in header)
    body = ''.join('<tr>' + ''.join(f'<td>{html.escape(str(value))}</td>' for value in row) + '</tr>' for row in rows)
    return f'<table id="{table_id or "table"}"><thead><tr>{headings}</tr></thead><tbody>{body}</tbody></table>'


def median_range(values, formatter=fmt_f):
    values = sorted(values)
    if not values:
        return '-'
    return formatter(values[len(values) // 2])


def status_cell(text, bad=False):
    return text


def box_traces(*args, **kwargs):
    return []


def sorted_bar(*args, **kwargs):
    return []


def vline(value):
    return {'type': 'line', 'x0': value, 'x1': value, 'y0': 0, 'y1': 1}


def alert_groups(groups):
    return [(name, values) for name, values in groups if values]


class Report:
    def __init__(self):
        self.sections = []

    def chart(self, *args, **kwargs):
        return 0

    def block(self, title, chart_ids, note=''):
        self.sections.append(f'<section><h2>{html.escape(title)}</h2><p>{html.escape(note)}</p></section>')

    def write(self, path, title):
        document = ('<!doctype html><html><head><meta charset="utf-8"><title>' + html.escape(title) +
                    '</title><style>body{font-family:sans-serif;max-width:1200px;margin:2rem auto}'
                    'table{border-collapse:collapse}th,td{border:1px solid #ccc;padding:.4rem}</style></head><body>' +
                    f'<h1>{html.escape(title)}</h1><p>Generated {datetime.now().isoformat(timespec="seconds")}</p>' +
                    ''.join(self.sections) + '</body></html>')
        with open(path, 'w', encoding='utf-8') as handle:
            handle.write(document)