import csv
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Extract iterations for GMG.')
    parser.add_argument('input_file', type=str, help='Input CSV file')
    parser.add_argument('smoothing', type=int, help='Smoothing value')
    parser.add_argument('shyness', type=float, help='Shyness value')
    parser.add_argument('threshold', type=float, help='Threshold value')
    args = parser.parse_args()

    writer = csv.writer(sys.stdout)
    writer.writerow(['p', 'refinements', 'iterations'])
    
    try:
        with open(args.input_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    row_shy = float(row['shyness'])
                    row_th = float(row['threshold'])
                    row_smooth = int(row['smoothing'])
                    
                    if (abs(row_shy - args.shyness) < 1e-9 and 
                        abs(row_th - args.threshold) < 1e-9 and 
                        row_smooth == args.smoothing):
                        
                        iterations = row.get('gmg_its')
                        if iterations and iterations.strip():
                            writer.writerow([row['p'], row['refinements'], iterations])
                except (ValueError, KeyError):
                    continue
    except FileNotFoundError:
        # Fail silently to avoid breaking LaTeX when invoked via -shell-escape
        sys.exit(0)

if __name__ == "__main__":
    main()
