#!/usr/bin/env python3
#file:settings.json
import csv
import sys
import os
import tempfile
import shutil
import json

# This script is intended to be run from the sbm_mg_patch_paper/ directory
IN_CSV = os.path.join(os.getcwd(), 'results', 'full_2D.csv')
# determine OUT_DIR from settings directive if present, otherwise default to results/tmp
default_tmp = 'results/tmp'
try:
    script_dir = os.path.dirname(__file__)
    try:
        with open(__file__, 'r') as sf:
            for _ in range(20):
                line = sf.readline()
                if not line:
                    break
                line = line.strip()
                if line.startswith('#file:'):
                    fname = line.split(':', 1)[1].strip()
                    settings_path = os.path.join(script_dir, fname)
                    if os.path.exists(settings_path):
                        with open(settings_path) as s2:
                            settings = json.load(s2)
                            default_tmp = settings.get('tmp_dir', default_tmp)
                    break
    except Exception:
        pass
except Exception:
    pass

if os.path.isabs(default_tmp):
    OUT_DIR = default_tmp
else:
    OUT_DIR = os.path.join(os.getcwd(), default_tmp)

OUT_CSV = os.path.join(OUT_DIR, 'full_2D.filtered.csv')

try:
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(IN_CSV, newline='') as f:
        r = csv.reader(f)
        try:
            header = next(r)
        except StopIteration:
            # create an empty filtered file with no rows
            open(OUT_CSV, 'w').close()
            sys.exit(0)
        if 'amg_time' not in header:
            # write header-only filtered file and exit
            with open(OUT_CSV, 'w', newline='') as out_f:
                csv.writer(out_f).writerow(header)
            sys.exit(0)
        idx = header.index('amg_time')
        out_rows = [header]
        for row in r:
            if len(row) <= idx:
                continue
            v = row[idx].strip()
            if v == '':
                continue
            try:
                float(v)
            except Exception:
                continue
            out_rows.append(row)
    # write to a temp file inside OUT_DIR then atomically move
    with tempfile.NamedTemporaryFile('w', delete=False, newline='', dir=OUT_DIR) as t:
        w = csv.writer(t)
        w.writerows(out_rows)
        tmp = t.name
    os.replace(tmp, OUT_CSV)
except Exception:
    # Fail silently to avoid breaking LaTeX when invoked via -shell-escape
    sys.exit(0)
