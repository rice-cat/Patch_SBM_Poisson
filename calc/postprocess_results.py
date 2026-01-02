#!/usr/bin/env python3
"""Postprocess results produced by runner.sh

Usage: postprocess_results.py RESULTS_DIR

This script scans immediate subdirectories of RESULTS_DIR for folders
matching the pattern produced by the runner, e.g.:

  2D_p3_dist0.01_001
  3D_p15_dist0.1_002

It extracts: dimension, FE degree, distortion, and run index, and writes
an `output.txt` CSV file under RESULTS_DIR with one row per found run.

Columns: folder, dim, degree, distortion, run_index, param_file, output_file
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from typing import Optional
import statistics


# New runner creates folders like: 2D_p1_ref5_shy3_smo3_th1.
PATTERN = re.compile(r'^(?P<dim>2D|3D|2d|3d)_p(?P<degree>\d+)_ref(?P<ref>\d+)(?:_shy(?P<shy>[^_]+))?(?:_smo(?P<smo>\d+))?(?:_th(?P<th>[^_]+))?$')


def parse_folder_name(name: str) -> Optional[dict]:
    m = PATTERN.match(name)
    if not m:
        return None
    return {
        'folder': name,
        'dim': m.group('dim'),
        'degree': int(m.group('degree')),
        'refinements': int(m.group('ref')),
        'shy': m.group('shy') or '',
        'smoothing': m.group('smo') or '',
        'threshold': m.group('th') or '',
    }


def find_runs(results_dir: str):
    entries = sorted(os.listdir(results_dir))
    for e in entries:
        path = os.path.join(results_dir, e)
        if not os.path.isdir(path):
            continue
        parsed = parse_folder_name(e)
        if not parsed:
            # skip non-matching folders
            continue
        # look for files inside folder
        param_file = os.path.join(path, 'used_parameters.prm')
        output_file = os.path.join(path, 'output.txt')
        if not os.path.exists(param_file):
            # try any .prm file in the folder
            for f in os.listdir(path):
                if f.lower().endswith('.prm'):
                    param_file = os.path.join(path, f)
                    break
        if not os.path.exists(output_file):
            raise RuntimeError(f'Missing output.txt in folder: {path}')

        parsed['param_file'] = param_file if os.path.exists(param_file) else ''
        parsed['output_file'] = output_file if os.path.exists(output_file) else ''

        # parse output file for CG iterations, and level-wise cells/DoFs
        parsed['iterations'] = ''
        parsed['cells_finest'] = ''
        parsed['dofs_finest'] = ''
        parsed['smoother_sweep_time'] = ''
        parsed['gmg_solve_time'] = ''
        parsed['gmg_iterations'] = ''
        parsed['amg_solve_time'] = ''
        parsed['amg_iterations'] = ''

        try:
            with open(output_file, 'r') as fh:
                text = fh.read()

            # iterations: look for "Solved in X iteration(s)" (allow singular/plural),
            # and optionally a trailing comma/period
            m_iter = re.search(r"Solved in\s*([0-9]+)\s+iteration(?:s)?\b[,\.]?", text)
            if m_iter:
                parsed['iterations'] = int(m_iter.group(1))
                if not parsed.get('gmg_iterations'):
                    parsed['gmg_iterations'] = int(m_iter.group(1))

            # Extract SUMMARY values
            m_sweep = re.search(r"SUMMARY: smoother_sweep_time\s+([0-9.eE+-]+)(?:\s*[a-zA-Z%]+)?", text)
            if m_sweep:
                parsed['smoother_sweep_time'] = float(m_sweep.group(1))

            m_gmg_time = re.search(r"SUMMARY: gmg_solve_time\s+([0-9.eE+-]+)(?:\s*[a-zA-Z%]+)?", text)
            if m_gmg_time:
                parsed['gmg_solve_time'] = float(m_gmg_time.group(1))
            else:
                # Fallback for older output format or missing SUMMARY
                # Look for "Solving system with multigrid" and then "Solved in X iterations"
                # We don't have a direct solve time without SUMMARY, but we might have it in a different line
                pass

            m_gmg_iter = re.search(r"SUMMARY: gmg_iterations\s+([0-9]+)", text)
            if m_gmg_iter:
                parsed['gmg_iterations'] = int(m_gmg_iter.group(1))

            m_amg_time = re.search(r"SUMMARY: amg_solve_time\s+([0-9.eE+-]+)(?:\s*[a-zA-Z%]+)?", text)
            if m_amg_time:
                parsed['amg_solve_time'] = float(m_amg_time.group(1))

            m_amg_iter = re.search(r"SUMMARY: amg_iterations\s+([0-9]+)", text)
            if m_amg_iter:
                parsed['amg_iterations'] = int(m_amg_iter.group(1))

            # Extract cell threshold if not already set from folder name
            if not parsed.get('threshold'):
                m_th_out = re.search(r"cell threshold\s+([0-9.]+)", text, re.I)
                if m_th_out:
                    parsed['threshold'] = m_th_out.group(1)
                else:
                    # Try to find it in the parameter file
                    if parsed.get('param_file'):
                        try:
                            with open(parsed['param_file'], 'r') as pf:
                                ptext = pf.read()
                            m_th_prm = re.search(r'^\s*set\s+cell_threshold\s*=\s*([^\s#]+)', ptext, re.I | re.M)
                            if m_th_prm:
                                parsed['threshold'] = m_th_prm.group(1)
                        except:
                            pass
                if not parsed.get('threshold'):
                    parsed['threshold'] = '1.'
                
            # Extract shyness if not already set from folder name
            if not parsed.get('shy'):
                m_shy = re.search(r"Shyness:\s+([0-9]+)", text)
                if m_shy:
                    parsed['shy'] = m_shy.group(1)
                elif parsed.get('param_file'):
                    try:
                        with open(parsed['param_file'], 'r') as pf:
                            ptext = pf.read()
                        m_shy_prm = re.search(r'^\s*set\s+shyness\s*=\s*([^\s#]+)', ptext, re.I | re.M)
                        if m_shy_prm:
                            parsed['shy'] = m_shy_prm.group(1)
                    except:
                        pass

            # Extract smoothing if not already set from folder name
            if not parsed.get('smoothing'):
                m_smo = re.search(r"smoothing steps:\s+([0-9]+)", text, re.I)
                if m_smo:
                    parsed['smoothing'] = m_smo.group(1)
                elif parsed.get('param_file'):
                    try:
                        with open(parsed['param_file'], 'r') as pf:
                            ptext = pf.read()
                        m_smo_prm = re.search(r'^\s*set\s+n_smoothing_steps\s*=\s*([^\s#]+)', ptext, re.I | re.M)
                        if m_smo_prm:
                            parsed['smoothing'] = m_smo_prm.group(1)
                    except:
                        pass

            # Find all "Level <n>: <num> cells" and "Level <n>: <num> DoFs"
            cells = {}
            dofs = {}
            for m in re.finditer(r"Level\s+(\d+):\s*([0-9]+)\s+cells", text):
                lvl = int(m.group(1)); val = int(m.group(2)); cells[lvl] = val
            for m in re.finditer(r"Level\s+(\d+):\s*([0-9]+)\s+DoFs", text):
                lvl = int(m.group(1)); val = int(m.group(2)); dofs[lvl] = val

            if cells:
                finest = max(cells.keys())
                parsed['cells_finest'] = cells.get(finest, '')
            if dofs:
                finest_d = max(dofs.keys())
                parsed['dofs_finest'] = dofs.get(finest_d, '')

        except Exception:
            # ignore parse errors and leave fields empty
            pass

        yield parsed

def write_csv(results_dir: str, rows):
    out_path = os.path.join(results_dir, 'output.txt')
    fieldnames = [
        'folder', 'dim', 'degree', 'refinements', 'shy', 'smoothing', 'threshold',
        'gmg_iterations', 'gmg_solve_time', 'smoother_sweep_time',
        'amg_iterations', 'amg_solve_time'
    ]
    # Collect all rows and convert fields
    formatted_rows = []
    for r in rows:
        row = {k: r.get(k, '') for k in fieldnames}
        dim_val = row['dim']
        if isinstance(dim_val, str) and dim_val.lower().startswith('2'):
            row['dim'] = '2'
        elif isinstance(dim_val, str) and dim_val.lower().startswith('3'):
            row['dim'] = '3'
        formatted_rows.append(row)
    # Compute max width for each column
    col_widths = {k: max(len(str(row[k])) for row in ([{k: k for k in fieldnames}] + formatted_rows)) for k in fieldnames}
    with open(out_path, 'w', newline='') as fh:
        # Write header
        fh.write('  '.join(f"{k:<{col_widths[k]}}" for k in fieldnames) + '\n')
        count = 0
        for row in formatted_rows:
            fh.write('  '.join(f"{str(row[k]):<{col_widths[k]}}" for k in fieldnames) + '\n')
            count += 1
    return out_path, count


def write_summary_csv(results_dir: str, rows):
    """Write a CSV summary (summary.csv) with one row per run.

    Columns: p, refinements, cells, dofs, iterations, ...
    """
    out_path = os.path.join(results_dir, 'summary.csv')
    fieldnames = [
        'p', 'refinements', 'shyness', 'smoothing', 'threshold', 'cells', 'dofs',
        'gmg_its', 'gmg_time', 'sweep_time', 'amg_its', 'amg_time'
    ]
    written = 0
    with open(out_path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(fieldnames)
        for r in rows:
            p = r.get('degree', '')
            ref = r.get('refinements', '')
            cells = r.get('cells_finest', '')
            dofs = r.get('dofs_finest', '')
            
            gmg_its = r.get('gmg_iterations', '')
            gmg_time = r.get('gmg_solve_time', '')
            sweep_time = r.get('smoother_sweep_time', '')
            amg_its = r.get('amg_iterations', '')
            amg_time = r.get('amg_solve_time', '')

            # Write rows even when GMG timing/iterations are missing (e.g., smoothing-only runs)

            # Try to read shyness, n_smoothing_steps, and cell_threshold from parameter file if available
            shyness = r.get('shy', '')
            smoothing = r.get('smoothing', '')
            threshold = r.get('threshold', '')
            
            writer.writerow([
                p, ref, shyness, smoothing, threshold, cells, dofs,
                gmg_its, gmg_time, sweep_time, amg_its, amg_time
            ])
            written += 1
    return out_path, written


def main(argv=None):
    ap = argparse.ArgumentParser(description='Postprocess results produced by runner.sh')
    ap.add_argument('results_dir', help='Path to the RESULTS directory produced by runner.sh')
    args = ap.parse_args(argv)

    results_dir = args.results_dir
    if not os.path.isdir(results_dir):
        print(f'ERROR: not a directory: {results_dir}', file=sys.stderr)
        return 2

    rows = list(find_runs(results_dir))
    if not rows:
        print('No matching run folders found under', results_dir)
        return 0

    out_path, count = write_csv(results_dir, rows)
    print(f'Wrote {count} rows to {out_path}')
    # write summary grouped by degree and distortion
    sum_path, sum_count = write_summary_csv(results_dir, rows)
    print(f'Wrote {sum_count} summary rows to {sum_path}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
