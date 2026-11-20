#!/usr/bin/env python3
"""Inspect a prospective PUBLIC release tree for private artifacts and identity leaks.
This is a supplementary check, not a substitute for human disclosure review.
"""
from pathlib import Path
import sys
root = Path(__file__).resolve().parents[1]
blocked_dirs = {'validation_private', 'results_private', 'raw', 'processed', '.Rproj.user'}
blocked_suffix = {'.xlsx', '.xls', '.rar', '.zip', '.fasta', '.bam', '.fastq'}
blocked_names = {'PROTEINAS_ORIGINALES.csv'}
problems=[]
for f in root.rglob('*'):
    if not f.is_file(): continue
    rel=f.relative_to(root)
    if any(part in blocked_dirs for part in rel.parts) or f.suffix.lower() in blocked_suffix or f.name in blocked_names or rel.as_posix() == 'data/peaks_manifest.csv':
        problems.append(str(rel))
if problems:
    print('FAIL: files not suitable for a public release:')
    print('\n'.join(problems))
    sys.exit(1)
print('PASS: no blocked private files found (human review still required).')
