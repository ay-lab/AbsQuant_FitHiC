import numpy as np
import pandas as pd
import subprocess as sp
import time
from multiprocessing import Pool
import cv2
import cooler
import argparse
import sys
import os

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Filter loops per chromosome')
    
    parser.add_argument('--per_replicate_files',
                       required=True,
                       nargs='+',
                       help='List of file paths to per-replicate .mcool files')
    
    parser.add_argument('--combined_replicate_file',
                       required=True,
                       help='File path to combined-replicate .mcool file')
    
    parser.add_argument('--P_s_curves_dir',
                       required=True,
                       help='Directory containing P(s) curves for per-replicate and combined replicate files')
    
    parser.add_argument('--fithic_combined_loop_file',
                       required=True,
                       help='File path to combined fithic loop call file (e.g., significances.txt.gz)')
    
    parser.add_argument('--per_replicate_names',
                       required=True,
                       nargs='+',
                       help='List of replicate names (e.g., npTh17-1 npTh17-2)')
    
    parser.add_argument('--combined_replicate_name',
                       required=True,
                       help='Name of the combined replicate (e.g., npTh17)')
    
    parser.add_argument('--resolution',
                       type=int,
                       default=10000,
                       help='Resolution in base pairs (default: 10000)')
    
    parser.add_argument('--nproc',
                       type=int,
                       default=30,
                       help='Number of processors to use (default: 30)')
    
    parser.add_argument('--chromosome',
                       required=True,
                       help='Chromosome to process (e.g., chr1)')
    
    parser.add_argument('--chunk_size',
                       type=int,
                       default=40,
                       help='Chunk size for multiprocessing (default: 40)')
    
    parser.add_argument('--output_dir',
                       required=True,
                       help='Output directory')
    
    parser.add_argument('--fdr_threshold',
                       type=float,
                       default=0.001,
                       help='FDR threshold for significance (default: 0.001)')
    
    parser.add_argument('--verbose',
                       action='store_true',
                       default=True,
                       help='Enable verbose output (default: True)')
    
    return parser.parse_args()

def acceptable_size_and_location(clr, chrom, left, right, local_region_size, min_size=32000):
    '''Check that loop is an appropriate size (>=min_size) and is not too close to end of chromosome.'''
    if right-left < min_size:
        return False
    if left-local_region_size < 0 or right+local_region_size > clr.chromsizes[chrom]:
        return False
    return True

def no_NaNs_near_center(clr, chrom, left, right, local_region_size, na_stripe_dist_to_center_px_cutoff=5):
    '''Check that there are no NaN stripes too close to the center.
    Args:
        clr: cooler object
        chrom: chromosome name
        left: left position of loop
        right: right position of loop
        local_region_size: size of local region in base pairs
        na_stripe_dist_to_center_px_cutoff: maximum distance from center to NaN stripe in pixels
    Returns:
        True if no NaN stripes too close to the center, False otherwise
    '''
    img = clr.matrix().fetch(f'{chrom}:{left-local_region_size}-{left+local_region_size}',
                            f'{chrom}:{right-local_region_size}-{right+local_region_size}').astype('float')
    
    length_of_img = img.shape[0]
    ver_na_stripe_indices = np.where(np.sum(np.isnan(img),0)==length_of_img)[0]
    hor_na_stripe_indices = np.where(np.sum(np.isnan(img),1)==length_of_img)[0]
    any_na_stripes = len(ver_na_stripe_indices)>0 or len(hor_na_stripe_indices)>0

    if any_na_stripes:
        middle_index = length_of_img//2
        ver_na_stripe_indices_from_middle = np.abs(ver_na_stripe_indices - middle_index)
        hor_na_stripe_indices_from_middle = np.abs(hor_na_stripe_indices - middle_index)
        if np.any(ver_na_stripe_indices_from_middle<=na_stripe_dist_to_center_px_cutoff) or \
           np.any(hor_na_stripe_indices_from_middle<=na_stripe_dist_to_center_px_cutoff):
            return False
    return True

def read_counts_per_pixel(clr, chrom, left, right, local_region_size):
    '''Calculate reads per pixel in local region.
    Args:
        clr: cooler object
        chrom: chromosome name
        left: left position of loop
        right: right position of loop
        local_region_size: size of local region in base pairs
    Returns:
        read_count: reads per pixel in local region
        num_pixels: number of pixels in local region
    '''
    img_unbalanced = clr.matrix(balance=False).fetch(
        f'{chrom}:{left-local_region_size}-{left+local_region_size}',
        f'{chrom}:{right-local_region_size}-{right+local_region_size}').astype('float')
    read_count = np.sum(img_unbalanced)
    num_pixels = np.size(img_unbalanced)
    return read_count/num_pixels

def global_maximum_dist_to_center(clr, chrom, left, right, P_s_data, s_px_matrix, local_region_size, x_px, y_px, gaussian_blur_sigma_px=2.5, ignore_diag_cutoff_px=5):
    '''
    Calculate the Euclidean distance (in pixels) of the global maximum to the center of the image (the location of the loop). The global maximum is calculated on the observed/expected matrix.
    Args:
        clr: cooler object
        chrom: chromosome name
        left: left position of loop
        right: right position of loop
        P_s_data: P(s) curves data
        s_px_matrix: s_px_matrix
        local_region_size: size of local region in base pairs
        x_px: x-coordinates of pixels
        y_px: y-coordinates of pixels
        gaussian_blur_sigma_px: sigma for Gaussian blur
        ignore_diag_cutoff_px: maximum distance from center to diagonal in pixels
    Returns:
        dist_to_brightest_pixel: Euclidean distance (in pixels) of the global maximum to the center of the image
    '''
    # get the image
    img = clr.matrix().fetch(f'{chrom}:{left-local_region_size}-{left+local_region_size}',f'{chrom}:{right-local_region_size}-{right+local_region_size}').astype('float')

    # get the expected global background image
    bg_img = P_s_data[s_px_matrix]

    # in the image and background image, make all pixels near diagonal NA
    img[s_px_matrix<=ignore_diag_cutoff_px] = np.nan
    bg_img[s_px_matrix<=ignore_diag_cutoff_px] = np.nan

    # divide the image by the expected global background
    img_over_bg = img/bg_img

    # resolve any NA values (only do this if no NaNs near center)
    img_over_bg_NAs_removed = np.nan_to_num(img_over_bg, nan=np.nanmedian(img))  # replace NA values with median value in the image
    
    # blur image
    ksize = int(np.ceil(3*gaussian_blur_sigma_px)//2*2+1)  # round up to next odd integer >= 3 sigma
    img_over_bg_blurred = cv2.GaussianBlur(img_over_bg_NAs_removed,ksize=(ksize,ksize),sigmaX=gaussian_blur_sigma_px)
    
    # find global maximum
    center_pixel_indices = np.array([i[0] for i in np.where(np.logical_and(x_px==0,y_px==0))])
    brightest_pixel_indices = np.array(np.unravel_index(np.nanargmax(img_over_bg_blurred), img_over_bg_blurred.shape))
    dist_to_brightest_pixel = np.linalg.norm(brightest_pixel_indices-center_pixel_indices)

    return dist_to_brightest_pixel

def get_image(clr, chrom, left, right, local_region_size, s_px_matrix, P_s_data=None, over_background=False, ignore_diag_cutoff_px=5):
    '''
    Get the image (for diagnostic purposes).
    Args:
        clr: cooler object
        chrom: chromosome name
        left: left position of loop
        right: right position of loop
        P_s_data: P(s) curves data
        over_background: whether to divide the image by the expected global background
        s_px_matrix: s_px_matrix
        local_region_size: size of local region in base pairs
        ignore_diag_cutoff_px: maximum distance from center to diagonal in pixels
    Returns:
        img: image (if over_background==False) or image over background (if over_background==True)
    '''
        
    # get the image
    img = clr.matrix().fetch(f'{chrom}:{left-local_region_size}-{left+local_region_size}',f'{chrom}:{right-local_region_size}-{right+local_region_size}').astype('float')
    
    # make all pixels near diagonal NA
    img[s_px_matrix<=ignore_diag_cutoff_px] = np.nan
        
    if over_background:
        # get the expected global background image
        bg_img = P_s_data[s_px_matrix]
        
        # make all pixels near diagonal NA
        bg_img[s_px_matrix<=ignore_diag_cutoff_px] = np.nan

        # divide the image by the expected global background
        img_over_bg = img/bg_img
        
        return img_over_bg
    else:
        return img

def load_coolers_and_P_s_curves(args):
    # Initialize coolers dictionary
    coolers = {}
    sample_names = []

    # Validate that number of names matches number of files
    if len(args.per_replicate_names) != len(args.per_replicate_files):
        print(f"Error: Number of per-replicate names ({len(args.per_replicate_names)}) does not match number of files ({len(args.per_replicate_files)})")
        sys.exit(1)

    # Load per-replicate coolers
    if args.verbose:
        print("Loading per-replicate cooler files...")
    for i, replicate_file in enumerate(args.per_replicate_files):
        sample_name = args.per_replicate_names[i]
        
        if args.verbose:
            print(f"  Loading: {replicate_file} (name: {sample_name})")
        
        if not os.path.exists(replicate_file):
            print(f"  ✗ Warning: Cooler file not found: {replicate_file}")
            continue
        
        try:
            cooler_path = f'{replicate_file}::/resolutions/{args.resolution}'
            coolers[sample_name] = cooler.Cooler(cooler_path)
            sample_names.append(sample_name)
            if args.verbose:
                print(f"  ✓ Successfully loaded: {sample_name}")
        except Exception as e:
            print(f"  ✗ Error loading cooler {replicate_file}: {e}")
            continue

    # Load combined-replicate cooler
    if args.verbose:
        print("")
        print("Loading combined replicate cooler file...")
    
    combined_sample_name = args.combined_replicate_name
    
    if args.verbose:
        print(f"  Loading: {args.combined_replicate_file}")
    
    if not os.path.exists(args.combined_replicate_file):
        print(f"  ✗ Warning: Combined cooler file not found: {args.combined_replicate_file}")
    else:
        try:
            combined_cooler_path = f'{args.combined_replicate_file}::/resolutions/{args.resolution}'
            coolers[combined_sample_name] = cooler.Cooler(combined_cooler_path)
            if args.verbose:
                print(f"  ✓ Successfully loaded: {combined_sample_name}")
        except Exception as e:
            print(f"  ✗ Error loading combined cooler {args.combined_replicate_file}: {e}")
    
    if len(coolers) == 0:
        print("Error: No cooler files were successfully loaded")
        sys.exit(1)

    # Load P(s) curves for each replicate
    if args.verbose:
        print("")
        print("Loading P(s) curves...")
    P_s_curves = {}
    for sample_name in coolers.keys():
        P_s_curve_path = f'{args.P_s_curves_dir}/{sample_name}.P_s_{args.resolution}bp.txt'
        if os.path.exists(P_s_curve_path):
            P_s_curves[sample_name] = np.loadtxt(P_s_curve_path)
            if args.verbose:
                print(f"  ✓ Loaded P(s) curve: {sample_name}")
        else:
            print(f"  ✗ Warning: P(s) curve not found: {P_s_curve_path}")
    
    return coolers, P_s_curves, sample_names, combined_sample_name


if __name__ == "__main__":
    args = parse_arguments()
    verbose = args.verbose
    
    if verbose:
        print("==================================================")
        print("Filtering Loops per Chromosome")
        print("==================================================")
        print(f"Per-replicate files: {len(args.per_replicate_files)} files")
        print(f"Combined replicate file: {args.combined_replicate_file}")
        print(f"P(s) curves directory: {args.P_s_curves_dir}")
        print(f"Fithic combined loop file: {args.fithic_combined_loop_file}")
        print(f"Chromosome: {args.chromosome}")
        print(f"Resolution: {args.resolution} bp")
        print(f"FDR threshold: {args.fdr_threshold}")
        print(f"Number of processors: {args.nproc}")
        print(f"Output directory: {args.output_dir}")
        print("")
    
    # Load coolers and P(s) curves
    coolers, P_s_curves, per_replicate_names, combined_replicate_name = load_coolers_and_P_s_curves(args)
    res = args.resolution
    
    # Parameters
    min_read_counts_per_pixel = 0.4
    local_region_size = res * 10
    a = local_region_size//res
    y_px, x_px = np.meshgrid(np.arange(-a, a+1),np.arange(-a, a+1))

    # load table of loops from fithic
    fithic_loops = None
    if os.path.exists(args.fithic_combined_loop_file):
        if verbose:
            print(f"Loading fithic loops from: {args.fithic_combined_loop_file}")
        fithic_loops = pd.read_csv(args.fithic_combined_loop_file, sep='\t')
        fithic_loops['start1'] = fithic_loops['fragmentMid1'].astype(int) - (res//2)
        fithic_loops['end1'] = fithic_loops['fragmentMid1'].astype(int) + (res//2)
        fithic_loops['start2'] = fithic_loops['fragmentMid2'].astype(int) - (res//2)
        fithic_loops['end2'] = fithic_loops['fragmentMid2'].astype(int) + (res//2)

        # filter for chromosome
        fithic_loops = fithic_loops[fithic_loops['chr1'] == args.chromosome]
        fithic_loops = fithic_loops[fithic_loops['chr2'] == args.chromosome]
        # filter for q-value
        fithic_loops = fithic_loops[fithic_loops['q-value'] < args.fdr_threshold]

        # table of loops
        loops_bedpe = pd.DataFrame(fithic_loops[['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']].values)
        if verbose:
            print(f"  ✓ Loaded {len(loops_bedpe)} loops from fithic file")
    else:
        print(f"Error: Fithic combined file not found: {args.fithic_combined_loop_file}")
        sys.exit(1)
    
    # Convert to internal format
    loops_df = pd.DataFrame(columns=['chr','left','right','size'])
    loops_df['chr'] = loops_bedpe[0]
    loops_df['left'] = (loops_bedpe[1] + loops_bedpe[2])//2
    loops_df['right'] = (loops_bedpe[4] + loops_bedpe[5])//2
    loops_df['size'] = loops_df['right'] - loops_df['left']
    
    # Initialize results dataframe
    column_names = ['size_loc_pass'] + \
                  [f'NaN_pass_R{i+1}' for i in range(len(per_replicate_names))] + \
                  [f'NaN_pass_all_merged'] + \
                  [f'read_counts_per_pixel_R{i+1}' for i in range(len(per_replicate_names))] + \
                  [f'global_max_dist_all_merged']
    loop_filtering_criteria_df = pd.DataFrame(index=loops_df.index, columns=column_names)

    def run_loop_filtering(k):
        # get loop
        loop = loops_df.loc[k]
        chrom = loop['chr']
        left = loop['left']
        right = loop['right']
        
        # adjust left and right to be in bin centers
        left = left//res*res + res//2
        right = right//res*res + res//2

        # calculate s_px_matrix
        loop_size_px = right//res-left//res
        s_px_matrix = loop_size_px+y_px-x_px  # genomic separation in units of res
        s_px_matrix[s_px_matrix<0] = 0  # don't allow negative values of s
        # Ensure s_px_matrix values don't exceed P(s) data bounds
        max_s_index = len(P_s_curves[combined_replicate_name]) - 1
        s_px_matrix[s_px_matrix > max_s_index] = max_s_index
        
        # Check size and location
        size_loc_pass = acceptable_size_and_location(coolers[per_replicate_names[0]], chrom, left, right, local_region_size)
        
        # Check NaN regions for each replicate
        nan_passes = [no_NaNs_near_center(coolers[rep], chrom, left, right, local_region_size) 
                    for rep in per_replicate_names]
        
        # Check NaN regions for combined replicate
        nan_pass_all_merged = no_NaNs_near_center(coolers[combined_replicate_name], chrom, left, right, local_region_size)
        
        if not (size_loc_pass and all(nan_passes)):
            return [size_loc_pass] + nan_passes + [nan_pass_all_merged] + [None] * len(per_replicate_names) + [None]
        
        # Check read counts for each replicate
        read_counts = [read_counts_per_pixel(coolers[rep], chrom, left, right, local_region_size) 
                    for rep in per_replicate_names]
        
        # Check global maximum distance to center
        global_max_dist_all_merged = global_maximum_dist_to_center(coolers[combined_replicate_name], chrom, left, right, P_s_curves[combined_replicate_name], s_px_matrix, local_region_size, x_px, y_px)

        
        return [size_loc_pass] + nan_passes + [nan_pass_all_merged] + read_counts + [global_max_dist_all_merged]
        
    # Set up multiprocessing
    chunk_size = args.chunk_size
    num_chunks = int(np.ceil(len(loops_df)/chunk_size))
    chunk_starts = np.arange(num_chunks)*chunk_size
    chunk_ends = (np.arange(num_chunks)+1)*chunk_size
    chunk_ends[-1] = len(loops_df)
    output_path = args.output_dir
    os.makedirs(output_path, exist_ok=True)
    
    # Run filtering
    for chunk_index in np.arange(num_chunks):
        start = time.time()
        with Pool(args.nproc) as p:
            indices_in_chunk = np.arange(chunk_starts[chunk_index], chunk_ends[chunk_index])
            results_in_chunk = p.map(run_loop_filtering, indices_in_chunk)
            loop_filtering_criteria_df.loc[indices_in_chunk] = results_in_chunk
        end = time.time()
        print(f'chunk {chunk_index+1} of {num_chunks} ({end-start:.2f} s)')
        if chunk_index % 100 == 0:
            loop_filtering_criteria_df.to_csv(f'{output_path}/loop_filtering_criteria_intermediate.fdr{args.fdr_threshold}.{args.chromosome}.csv', sep='\t', index=False)
    loop_filtering_criteria_df.to_csv(f'{output_path}/loop_filtering_criteria.fdr{args.fdr_threshold}.{args.chromosome}.csv', sep='\t', index=False)
       
    # Apply filters and save filtered loops
    if verbose:
        print("")
        print("Applying filters and saving filtered loops...")
    filter_conditions = [loop_filtering_criteria_df['size_loc_pass']]
    for i in range(len(per_replicate_names)):
        filter_conditions.append(loop_filtering_criteria_df[f'NaN_pass_R{i+1}'])
        filter_conditions.append(loop_filtering_criteria_df[f'read_counts_per_pixel_R{i+1}'] > min_read_counts_per_pixel)
    filter_conditions.append(loop_filtering_criteria_df['NaN_pass_all_merged'])
    filter_conditions.append(loop_filtering_criteria_df['global_max_dist_all_merged'] <= 2.5)
    filtered_indices = loops_df.index[np.all(filter_conditions, axis=0)]
    filtered_chrom_df = loops_df.loc[filtered_indices]
    
    output_file = f'{output_path}/filtered_loops.fdr{args.fdr_threshold}.{args.chromosome}.txt'
    filtered_chrom_df.to_csv(output_file, sep='\t', index=False)
    
    if verbose:
        print(f"  ✓ Saved {len(filtered_chrom_df)} filtered loops to: {output_file}")
        print("")
        print("==================================================")
        print("Analysis complete!")
        print("==================================================")
