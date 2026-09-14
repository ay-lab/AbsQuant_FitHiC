import argparse
import os
import sys

import pandas as pd

# Columns that define a loop. Both the concatenated table from 1.3 and a fithic
# significances table carry these names, which is what makes the join possible.
JOIN_COLS = ['chr1', 'fragmentMid1', 'chr2', 'fragmentMid2']


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Intersect a filtered loop set with per-replicate fithic calls')

    parser.add_argument('--combined_loops_file',
                        required=True,
                        help='Concatenated loop table from 1.3 '
                             '(<subset>.coords.fdr<fdr>.txt)')

    parser.add_argument('--replicate_files',
                        required=True,
                        nargs='+',
                        help='One or more fithic significances files to intersect '
                             'against. Plain text or .gz; both are read directly.')

    parser.add_argument('--output_dir',
                        required=True,
                        help='Output directory')

    parser.add_argument('--fdr_threshold',
                        type=float,
                        default=0.001,
                        help='FDR threshold, used in the output filename and, with '
                             '--apply_fdr, to filter the replicate calls '
                             '(default: 0.001)')

    parser.add_argument('--apply_fdr',
                        action='store_true',
                        help='Also drop replicate rows with q-value >= --fdr_threshold. '
                             'Off by default: the loop set from 1.3 is already '
                             'FDR-filtered, and the replicate table is being asked '
                             '"is this loop present here", not "is it significant here".')

    parser.add_argument('--verbose',
                        action='store_true',
                        default=True,
                        help='Enable verbose output (default: True)')

    return parser.parse_args()


def output_path_for(replicate_file, fdr_threshold, output_dir):
    """Build the output filename from the replicate filename.

    <name>.<rest>.significances.txt[.gz]
        -> <output_dir>/<name>.<rest>.fdr<fdr>.significances.txt

    Matches the naming the R version produced, so downstream steps that glob for
    *.fdr<fdr>.significances.txt keep working.
    """
    base = os.path.basename(replicate_file)
    if base.endswith('.gz'):
        base = base[:-len('.gz')]

    parts = base.split('.')
    replicate_name = parts[0]
    extension = '.'.join(parts[1:])

    if 'significances.txt' in extension:
        extension = extension.replace(
            'significances.txt', f'fdr{fdr_threshold}.significances.txt', 1)
    else:
        # Not a standard fithic filename; keep it unambiguous rather than
        # silently writing over something.
        extension = f'{extension}.fdr{fdr_threshold}.txt' if extension \
            else f'fdr{fdr_threshold}.txt'

    return replicate_name, os.path.join(output_dir, f'{replicate_name}.{extension}')


def load_combined_loops(path, verbose):
    """Read the 1.3 output, keeping only the four coordinate columns."""
    if not os.path.exists(path):
        print(f"Error: Combined loops file not found: {path}")
        sys.exit(1)

    df = pd.read_csv(path, sep='\t')

    missing = [c for c in JOIN_COLS if c not in df.columns]
    if missing:
        print(f"Error: {path} is missing column(s): {', '.join(missing)}")
        print(f"  Found: {', '.join(df.columns)}")
        print("  Expected the output of 1.3_concat_loops.sh.")
        sys.exit(1)

    df = df[JOIN_COLS].copy()
    df['fragmentMid1'] = df['fragmentMid1'].astype('int64')
    df['fragmentMid2'] = df['fragmentMid2'].astype('int64')

    n_before = len(df)
    df = df.drop_duplicates()
    if verbose:
        print(f"  ✓ Loaded {n_before} loops"
              + (f" ({n_before - len(df)} duplicate coordinates dropped)"
                 if len(df) != n_before else ""))
    return df


if __name__ == "__main__":
    args = parse_arguments()
    verbose = args.verbose

    if verbose:
        print("==================================================")
        print("Intersecting loops with per-replicate fithic calls")
        print("==================================================")
        print(f"Combined loops file: {args.combined_loops_file}")
        print(f"Replicate files: {len(args.replicate_files)}")
        print(f"FDR threshold: {args.fdr_threshold}")
        print(f"Apply FDR to replicates: {args.apply_fdr}")
        print(f"Output directory: {args.output_dir}")
        print("")

    os.makedirs(args.output_dir, exist_ok=True)

    if verbose:
        print(f"Reading combined loops file: {args.combined_loops_file}")
    combined_df = load_combined_loops(args.combined_loops_file, verbose)
    if verbose:
        print("")

    n_ok = 0
    for replicate_file in args.replicate_files:
        if verbose:
            print(f"Reading replicate file: {replicate_file}")

        if not os.path.exists(replicate_file):
            print(f"  ✗ Warning: Replicate file not found: {replicate_file}")
            continue

        # pandas reads .gz transparently, so no decompression step is needed.
        try:
            replicate_df = pd.read_csv(replicate_file, sep='\t')
        except Exception as e:
            print(f"  ✗ Error reading {replicate_file}: {e}")
            continue

        missing = [c for c in JOIN_COLS if c not in replicate_df.columns]
        if missing:
            print(f"  ✗ Warning: missing column(s) {', '.join(missing)} - skipping")
            continue

        replicate_df['fragmentMid1'] = replicate_df['fragmentMid1'].astype('int64')
        replicate_df['fragmentMid2'] = replicate_df['fragmentMid2'].astype('int64')
        n_replicate = len(replicate_df)

        if args.apply_fdr:
            if 'q-value' in replicate_df.columns:
                replicate_df = replicate_df[replicate_df['q-value'] < args.fdr_threshold]
                if verbose:
                    print(f"  q-value < {args.fdr_threshold}: "
                          f"{n_replicate} -> {len(replicate_df)} rows")
            else:
                print("  ✗ Warning: --apply_fdr given but no 'q-value' column; "
                      "not filtering")

        # Inner join: keep only the loops present in both tables, carrying the
        # replicate's own statistics through.
        filtered_df = replicate_df.merge(combined_df, on=JOIN_COLS, how='inner')
        # Put the coordinate columns first, as the data.table version did.
        filtered_df = filtered_df[JOIN_COLS +
                                  [c for c in filtered_df.columns if c not in JOIN_COLS]]

        replicate_name, output_file = output_path_for(
            replicate_file, args.fdr_threshold, args.output_dir)
        filtered_df.to_csv(output_file, sep='\t', index=False)

        pct = (100.0 * len(filtered_df) / len(combined_df)) if len(combined_df) else 0.0
        if verbose:
            print(f"  ✓ {replicate_name}: {len(filtered_df)} of {len(combined_df)} "
                  f"loops recovered ({pct:.1f}%), from {n_replicate} replicate calls")
            print(f"  ✓ Saved: {output_file}")
            print("")
        n_ok += 1

    if verbose:
        print("==================================================")
        print(f"Intersected {n_ok} of {len(args.replicate_files)} replicate files")
        print("==================================================")

    if n_ok == 0:
        sys.exit(1)
