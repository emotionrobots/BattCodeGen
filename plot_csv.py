#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Plot selected CSV columns (Y) against the first column (X)."
    )
    parser.add_argument("csvfile", help="Path to the CSV file")
    parser.add_argument(
        "columns",
        nargs="+",
        help="Column indices or names to plot against the first column (X)"
    )
    parser.add_argument("--delimiter", default=",", help="CSV delimiter (default=comma)")
    parser.add_argument("--header", type=int, default=0,
                        help="Row index for header (default=0, use None for no header)")
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.csvfile, delimiter=args.delimiter, header=args.header)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        sys.exit(1)

    # First column is X
    x = df.iloc[:, 0]
    x_label = df.columns[0] if args.header is not None else "X"

    plt.figure(figsize=(10, 6))
    for col in args.columns:
        try:
            # Try interpreting as integer index first
            col_idx = int(col)
            y = df.iloc[:, col_idx]
            label = df.columns[col_idx] if args.header is not None else f"col{col_idx}"
        except ValueError:
            # Otherwise treat as column name
            if col not in df.columns:
                print(f"Column '{col}' not found.")
                continue
            y = df[col]
            label = col

        plt.plot(x, y, label=label)

    plt.xlabel(x_label)
    plt.ylabel("Value")
    plt.title("CSV Plotter")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()

