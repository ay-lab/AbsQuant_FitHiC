import cooler
import argparse
import sys
import os

def open_cooler(path, resolution, verbose=False):
    """Open a cooler at `resolution`, accepting either .mcool or .cool.

    A multi-resolution .mcool holds several matrices and needs the
    "::/resolutions/<binsize>" URI suffix. A single-resolution .cool is opened
    directly and already has a fixed binsize. Each form fails on the other, so
    the file is inspected rather than assumed from its extension.

    A path that already contains "::" is passed through untouched, so an
    explicit URI still works.

    Args:
        path: path to a .cool or .mcool file, or a full cooler URI
        resolution: requested bin size in bp
        verbose: whether to note single-resolution files
    Returns:
        cooler.Cooler
    Raises:
        ValueError: if a .cool's own binsize is not `resolution`. Silently
            analysing at a different resolution would corrupt the P(s) indexing,
            which is in units of bins.
    """
    if '::' in path:
        return cooler.Cooler(path)

    try:
        multires = cooler.fileops.is_multires_file(path)
    except Exception:
        # Older cooler, or an unreadable header: fall back to the extension.
        multires = path.endswith('.mcool')

    if multires:
        return cooler.Cooler(f'{path}::/resolutions/{resolution}')

    clr = cooler.Cooler(path)
    if clr.binsize != resolution:
        raise ValueError(
            f"{path} is a single-resolution cooler with binsize {clr.binsize}, "
            f"but the requested resolution is {resolution}. Either pass "
            f"--resolution {clr.binsize}, or use an .mcool containing {resolution}.")
    if verbose:
        print(f"    (single-resolution .cool, binsize {clr.binsize})")
    return clr
# end def

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Calculate P(s) curves for Hi-C data')
    
    # Add arguments
    parser.add_argument('--per_replicate_files',
                       required=True,
                       nargs='+',
                       help='List of file paths to per-replicate .cool or .mcool files')
    
    parser.add_argument('--combined_replicate_file',
                       required=True,
                       help='File path to the combined-replicate .cool or .mcool file')
    
    # looptools.py ships in this repo, next to this script, so default to the
    # script's own directory. The previous default pointed at a separate
    # AbsLoopQuant_analysis_code checkout, which made a clone of this repo fail
    # at `import looptools` anywhere but the machine it was written on.
    parser.add_argument('--looptools_path',
                       default=os.path.dirname(os.path.abspath(__file__)),
                       help='Path to the directory containing looptools.py '
                            '(default: the directory holding this script)')
    
    parser.add_argument('--resolution',
                       type=int,
                       default=10000,
                       help='Resolution in base pairs (default: 10000)')
    
    parser.add_argument('--nproc',
                       type=int,
                       default=30,
                       help='Number of processors to use (default: 30)')
    
    parser.add_argument('--output_dir',
                       required=True,
                       help='Output directory')
    
    parser.add_argument('--verbose',
                       action='store_true',
                       default=True,
                       help='Enable verbose output (default: True)')
    
    return parser.parse_args()

def main():
    # Parse arguments
    args = parse_arguments()
    verbose = args.verbose
    
    if verbose:
        print("==================================================")
        print("Calculating P(s) Curves for Hi-C Data")
        print("==================================================")
        print(f"Per-replicate files: {len(args.per_replicate_files)} files")
        print(f"Combined replicate file: {args.combined_replicate_file}")
        print(f"Looptools path: {args.looptools_path}")
        print(f"Resolution: {args.resolution} bp")
        print(f"Number of processors: {args.nproc}")
        print(f"Output directory: {args.output_dir}")
        print("")
    
    # Import looptools
    if verbose:
        print(f"Importing looptools from: {args.looptools_path}")
    sys.path.insert(1, args.looptools_path)
    import looptools
    
    # Initialize coolers dictionary
    coolers = {}
    
    # Load per-replicate coolers
    if verbose:
        print("")
        print("Loading per-replicate cooler files...")
    for replicate_file in args.per_replicate_files:
        # Extract sample name from file path
        sample_name = os.path.splitext(os.path.basename(replicate_file))[0]
        
        if verbose:
            print(f"  Loading: {replicate_file}")
        
        if not os.path.exists(replicate_file):
            print(f"  ✗ Warning: Cooler file not found: {replicate_file}")
            continue
        
        try:
            coolers[sample_name] = open_cooler(replicate_file, args.resolution, args.verbose)
            if verbose:
                print(f"  ✓ Successfully loaded: {sample_name}")
        except Exception as e:
            print(f"  ✗ Error loading cooler {replicate_file}: {e}")
            continue
    
    # Load combined replicate cooler
    if verbose:
        print("")
        print("Loading combined replicate cooler file...")
    
    # Extract sample name from combined file path
    combined_sample_name = os.path.splitext(os.path.basename(args.combined_replicate_file))[0]
    
    if verbose:
        print(f"  Loading: {args.combined_replicate_file}")
    
    if not os.path.exists(args.combined_replicate_file):
        print(f"  ✗ Warning: Combined cooler file not found: {args.combined_replicate_file}")
    else:
        try:
            coolers[combined_sample_name] = open_cooler(
                args.combined_replicate_file, args.resolution, args.verbose)
            if verbose:
                print(f"  ✓ Successfully loaded: {combined_sample_name}")
        except Exception as e:
            print(f"  ✗ Error loading combined cooler {args.combined_replicate_file}: {e}")
    
    if len(coolers) == 0:
        print("Error: No cooler files were successfully loaded")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Calculate P(s) curves
    if verbose:
        print("")
        print("Calculating P(s) curves...")
        print("")
    
    for sample_name, clr in coolers.items():
        output_path = f'{args.output_dir}/{sample_name}.P_s_{args.resolution}bp.txt'
        
        if verbose:
            print(f"Processing: {sample_name}")
        
        if os.path.exists(output_path):
            if verbose:
                print(f"  P(s) curve already exists: {output_path}")
            continue
        
        try:
            if verbose:
                print(f"  Calculating P(s) curve...")
            looptools.calculate_and_save_avg_Ps_curve(clr, nproc=args.nproc, output_filename=output_path)
            if verbose:
                print(f"  ✓ Saved P(s) curve to: {output_path}")
        except Exception as e:
            print(f"  ✗ Error calculating P(s) curve for {sample_name}: {e}")
            continue
    
    if verbose:
        print("")
        print("==================================================")
        print("Analysis complete!")
        print("==================================================")

if __name__ == "__main__":
    main()