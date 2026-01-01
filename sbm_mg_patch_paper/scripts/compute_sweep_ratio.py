import csv
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Compute ratio of timing for 2 sweeps vs 1 sweep.')
    parser.add_argument('input_file', type=str, help='Input CSV file')
    parser.add_argument('shyness', type=float, nargs='?', default=3.0,
                        help='Shyness value (optional). Defaults to 3')
    parser.add_argument('threshold', type=float, nargs='?', default=1.0,
                        help='Threshold value (optional). Defaults to 1.0')
    args = parser.parse_args()

    data = []
    try:
        with open(args.input_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    # Convert values to float for comparison
                    row_shy = float(row['shyness'])
                    row_th = float(row['threshold'])
                    shy_match = (args.shyness is None) or (abs(row_shy - args.shyness) < 1e-9)
                    th_match = (args.threshold is None) or (abs(row_th - args.threshold) < 1e-9)
                    if shy_match and th_match:
                        data.append({
                            'p': int(row['p']),
                            'refinements': int(row['refinements']),
                            'smoothing': int(row['smoothing']),
                            'gmg_its': int(row['gmg_its']),
                            'gmg_time': float(row['gmg_time'])
                        })
                except (ValueError, KeyError):
                    continue
    except FileNotFoundError:
        print(f"File not found: {args.input_file}", file=sys.stderr)
        return

    if not data:
        print(f"No data found for shyness={args.shyness} and threshold={args.threshold}", file=sys.stderr)
        return

    # Group by (p, refinements)
    results = {}
    for entry in data:
        key = (entry['p'], entry['refinements'])
        if key not in results:
            results[key] = {}
        results[key][entry['smoothing']] = entry['gmg_time'] / entry['gmg_its']

    # Output CSV
    writer = csv.writer(sys.stdout)
    writer.writerow(['p', 'refinements', 'ratio'])
    
    for (p, ref), smooths in sorted(results.items()):
        if 1 in smooths and 2 in smooths:
            ratio = smooths[2] / smooths[1]
            writer.writerow([p, ref, ratio])

if __name__ == "__main__":
    main()
