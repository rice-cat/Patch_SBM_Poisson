import csv
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Compute throughput (DOF/s) for GMG and AMG.')
    parser.add_argument('input_file', type=str, help='Input CSV file')
    parser.add_argument('smoothing', type=int, help='Smoothing value')
    parser.add_argument('shyness', type=float, help='Shyness value')
    parser.add_argument('threshold', type=float, help='Threshold value')
    args = parser.parse_args()

    writer = csv.writer(sys.stdout)
    writer.writerow(['p', 'refinements', 'gmg_throughput', 'amg_throughput'])
    
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
                        
                        dofs = float(row['dofs'])
                        
                        gmg_throughput = ""
                        gmg_time = row.get('gmg_time')
                        if gmg_time and gmg_time.strip():
                            gmg_throughput = dofs / float(gmg_time)
                        else:
                            # Fallback to sweep_time * gmg_its
                            sweep_time = row.get('sweep_time')
                            gmg_its = row.get('gmg_its')
                            if sweep_time and sweep_time.strip() and gmg_its and gmg_its.strip():
                                gmg_throughput = dofs / (float(sweep_time) * float(gmg_its))
                                
                        amg_throughput = ""
                        amg_time = row.get('amg_time')
                        if amg_time and amg_time.strip():
                            amg_throughput = dofs / float(amg_time)
                            
                        if gmg_throughput != "" or amg_throughput != "":
                            writer.writerow([row['p'], row['refinements'], gmg_throughput, amg_throughput])
                except (ValueError, KeyError):
                    continue
    except FileNotFoundError:
        print(f"File not found: {args.input_file}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
