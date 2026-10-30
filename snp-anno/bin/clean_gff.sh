#!/usr/bin/env bash
set -euo pipefail
in="$1"
case "$in" in *.gz) reader=(zcat "$in") ;; *) reader=(cat "$in") ;; esac
"${reader[@]}" | tr -d '\r' | awk -F'\t' 'BEGIN { OFS = "\t" }
    /^#/ || NF < 9 { print; next }
    { n = split($9, attrs, ";"); out = ""; for (i = 1; i <= n; i++)
      if (attrs[i] !~ /^transcript_id=/) out = out (out == "" ? "" : ";") attrs[i]; $9 = out; print }'