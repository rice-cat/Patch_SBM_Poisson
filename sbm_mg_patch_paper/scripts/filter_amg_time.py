#!/usr/bin/env python3
"""Filter CSV rows keeping only those with a numeric `amg_time`.

Usage: filter_amg_time.py <input_csv> [output_dir]

Writes <output_dir>/<basename>.filtered.csv (defaults to results/tmp).
Exits silently (non-zero exit) on errors to avoid breaking LaTeX builds.
"""
import argparse
import csv
import json
import os
import sys
import tempfile


def main():
    p = argparse.ArgumentParser(description="Filter CSV by numeric amg_time column")
    p.add_argument('input_csv', help='Path to input CSV file')
    p.add_argument('output_dir', nargs='?', default=None, help='Directory for filtered output (overrides settings)')
    args = p.parse_args()

    in_path = args.input_csv
    out_dir = args.output_dir
    # load settings (script-local) for default tmp dir
    try:
        script_dir = os.path.dirname(__file__)
        settings_path = os.path.join(script_dir, 'settings.json')
        if os.path.exists(settings_path):
            with open(settings_path) as sf:
                settings = json.load(sf)
                default_tmp = settings.get('tmp_dir', 'results/tmp')
        else:
            default_tmp = 'results/tmp'
    except Exception:
        default_tmp = 'results/tmp'

    if out_dir is None:
        out_dir = default_tmp

    try:
        os.makedirs(out_dir, exist_ok=True)
        with open(in_path, newline='') as f:
            reader = csv.reader(f)
            try:
                header = next(reader)
            except StopIteration:
                # create empty output and exit
                open(os.path.join(out_dir, os.path.basename(in_path) + '.filtered.csv'), 'w').close()
                return
            if 'amg_time' not in header:
                # write header-only filtered file
                out_path = os.path.join(out_dir, os.path.basename(in_path) + '.filtered.csv')
                with open(out_path, 'w', newline='') as out_f:
                    csv.writer(out_f).writerow(header)
                return
            idx = header.index('amg_time')
            out_rows = [header]
            for row in reader:
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

        # use basename without extension for output file: e.g. full_2D.filtered.csv
        base = os.path.splitext(os.path.basename(in_path))[0]
        out_path = os.path.join(out_dir, base + '.filtered.csv')
        # write atomically in the output dir
        with tempfile.NamedTemporaryFile('w', delete=False, newline='', dir=out_dir) as t:
            writer = csv.writer(t)
            writer.writerows(out_rows)
            tmp = t.name
        os.replace(tmp, out_path)
    except Exception:
        # Fail silently to avoid breaking LaTeX when invoked via -shell-escape
        try:
            sys.exit(0)
        except SystemExit:
            return


if __name__ == '__main__':
    main()
