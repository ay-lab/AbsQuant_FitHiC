# AbsLoopQuant_TB - Absolute Loop Quantification Toolkit

This package provides tools for calculating P(s) curves from Hi-C data and filtering loops based on quantitative criteria. It implements the AbsLoopQuant approach for identifying high-confidence chromatin loops.

## Overview

The AbsLoopQuant workflow consists of two main steps:

1. **Calculate P(s) curves** - Compute probability of contact as a function of genomic distance from Hi-C cooler files
2. **Filter loops** - Filter loop calls based on size, location, NaN regions, read counts, and global maximum distance criteria

## Files

- **`1.1_calculate_P_s_curves_general.py`** - Calculate P(s) curves from Hi-C cooler files
- **`1.1_calculate_P_s_curves_general.sh`** - Bash script to generate SLURM jobs for P(s) curve calculation
- **`1.2_filter_loops.py`** - Filter loops per chromosome using quantitative criteria
- **`1.2_filter_loops.sh`** - Bash script to generate SLURM jobs for loop filtering
- **`looptools.py`** - Helper module with loop analysis utilities
- **`create_absloopquantTB_env.sh`** - Builds the `absloopquantTB` mamba environment and verifies its imports
- **`absloopquantTB_env.yml`** - Portable environment spec (any platform)
- **`absloopquantTB_env.lock.yml`** - Exact build-pinned export of the validated env (linux-64 only)
- **`archive/export_absloopquantTB_env.sh`** - Regenerates the lock file from a built env

## Setup

### 1. Create the environment

```bash
bash create_absloopquantTB_env.sh
```

This solves `absloopquantTB_env.yml`, then verifies that `numpy`, `pandas`,
`scipy`, `cooler`, `cooltools`, `cv2`, `matplotlib` and `looptools` all import -
a successful solve alone is not proof the env works.

Two specs are provided:

| File | Use it when |
|---|---|
| `absloopquantTB_env.yml` | **Default.** Portable; pins only what is load-bearing, solves on any platform. |
| `absloopquantTB_env.lock.yml` | You need the exact env validated on the LJI cluster. 349 build-pinned packages, **linux-64 only**. `bash create_absloopquantTB_env.sh --lock` |

Two pins in the portable spec are not cosmetic:

- **`numpy<2`** - `cooler`, `cooltools`, `py-opencv` and `numba` ship C
  extensions built against the numpy 1.x ABI, and numpy 2 breaks them at import.
- **`matplotlib-base<3.9`** - `cooltools.lib.plotting` imports `register_cmap`,
  removed in matplotlib 3.9, and its fallback targets a submodule that has never
  existed. `looptools` imports that module at load time, so **step 1 fails at
  import** without this pin.

Options:

```bash
ENV_NAME=absloopquantTB2 bash create_absloopquantTB_env.sh      # side-by-side build
ENV_PREFIX=/path/to/envs/absloopquantTB bash create_absloopquantTB_env.sh   # prefix env
```

Prefer `ENV_PREFIX` when home is small, quota'd, or unreliable.

### 2. Activate the Environment

```bash
mamba activate absloopquantTB
```

## Workflow

### Step 1: Calculate P(s) Curves

P(s) curves represent the probability of contact as a function of genomic distance. These are used as background models for loop filtering.

#### Generate SLURM Scripts

```bash
bash 1.1_calculate_P_s_curves_general.sh
```

This will:
- Find all per-replicate and combined `.mcool` files for each subset
- Generate SLURM scripts in `qshs/{date}_calculate_P_s_curves/`
- Create one script per subset

#### Run Python Script Directly

```bash
python3 1.1_calculate_P_s_curves_general.py \
    --per_replicate_files rep1.mcool rep2.mcool \
    --combined_replicate_file combined.mcool \
    --per_replicate_names rep1 rep2 \
    --combined_replicate_name combined \
    --output_dir /path/to/P_s_curves \
    --resolution 10000 \
    --verbose
```

#### Arguments

**Required:**
- `--per_replicate_files`: List of per-replicate `.mcool` file paths
- `--combined_replicate_file`: Path to combined replicate `.mcool` file
- `--per_replicate_names`: List of per-replicate names
- `--combined_replicate_name`: Name of combined replicate
- `--output_dir`: Output directory for P(s) curves

**Optional:**
- `--looptools_path`: Directory containing `looptools.py` (default: the directory holding the script, i.e. this repo)
- `--resolution`: Resolution in base pairs (default: 10000)
- `--nproc`: Number of processors (default: 30)
- `--verbose`: Enable verbose output (default: True)

Step 1 is resumable but not self-correcting: it skips any sample whose output
file already exists. After a job is killed mid-write, delete the partial
`.P_s_*bp.txt` before rerunning or the truncated file is kept silently.

#### Output

P(s) curve files are saved as:
```
{output_dir}/{sample_name}.P_s_{resolution}bp.txt
```

This is exactly the path step 2 reads as
`{P_s_curves_dir}/{sample_name}.P_s_{resolution}bp.txt`, so point
`--P_s_curves_dir` in step 2 at step 1's `--output_dir`.

### Step 2: Filter Loops

Filter loop calls from fithic based on quantitative criteria.

#### Generate SLURM Scripts

```bash
bash 1.2_filter_loops.sh
```

This will:
- Find all per-replicate and combined `.mcool` files
- Find fithic loop files
- Generate SLURM scripts in `qshs/{date}_filter_loops_per_chr_fdr{fdr}/`
- Create one script per chromosome per subset

#### Run Python Script Directly

```bash
python3 1.2_filter_loops.py \
    --per_replicate_files rep1.mcool rep2.mcool \
    --combined_replicate_file combined.mcool \
    --per_replicate_names rep1 rep2 \
    --combined_replicate_name combined \
    --P_s_curves_dir /path/to/P_s_curves \
    --fithic_combined_loop_file /path/to/significances.txt.gz \
    --chromosome chr1 \
    --output_dir /path/to/output \
    --resolution 10000 \
    --fdr_threshold 0.01 \
    --nproc 30 \
    --verbose
```

#### Arguments

**Required:**
- `--per_replicate_files`: List of per-replicate `.mcool` file paths
- `--combined_replicate_file`: Path to combined replicate `.mcool` file
- `--per_replicate_names`: List of per-replicate names
- `--combined_replicate_name`: Name of combined replicate
- `--P_s_curves_dir`: Directory containing P(s) curves
- `--fithic_combined_loop_file`: Path to fithic combined loop file (e.g., `significances.txt.gz`)
- `--chromosome`: Chromosome to process (e.g., `chr1`)
- `--output_dir`: Output directory for filtered loops

**Optional:**
- `--resolution`: Resolution in base pairs (default: 10000)
- `--fdr_threshold`: FDR threshold for significance (default: 0.001)
- `--nproc`: Number of processors (default: 30)
- `--chunk_size`: Chunk size for multiprocessing (default: 40)
- `--verbose`: Enable verbose output (default: True)

#### Filtering Criteria

Loops must pass all of the following criteria:

1. **Size and location**: Loop size ≥ 32,000 bp and not too close to chromosome ends
2. **NaN regions**: No NaN stripes too close to the center in any replicate
3. **Read counts**: Read counts per pixel ≥ 0.4 in each per-replicate
4. **Global maximum distance**: Distance of global maximum to center ≤ 2.5 pixels (in combined replicate)

#### Output

- **Filtering criteria file**: `loop_filtering_criteria.fdr{fdr}.{chromosome}.csv` - Contains all filtering criteria for each loop
- **Filtered loops file**: `filtered_loops.fdr{fdr}.{chromosome}.txt` - Final filtered loops passing all criteria

## Configuring the job generators

`1.1_calculate_P_s_curves_general.sh` and `1.2_filter_loops.sh` resolve this
repo from their own location, so a clone works from any path. The data project
is site-specific and is set by environment variable rather than by editing the
files:

```bash
BASE_DIR=/mnt/.../projects/<your-project> bash 1.2_filter_loops.sh
WORKING_DIR=/path/to/AbsLoopQuant_TB bash 1.2_filter_loops.sh   # code elsewhere
```

The `subsets`, `resolution`, `nproc` and `fdrThreshold` values near the top of
each generator still need editing per project.

## Known upstream bug, fixed here

`1.2_filter_loops.py` carries a fix that is **not** in the original
AbsLoopQuant code. `acceptable_size_and_location` correctly rejects a loop whose
±`local_region_size` window runs past the end of a chromosome, but upstream
only acts on that verdict *after* calling `no_NaNs_near_center`, which fetches
the out-of-bounds region first:

```python
size_loc_pass = acceptable_size_and_location(...)      # correctly False
nan_passes = [no_NaNs_near_center(...) for rep in ...] # runs anyway -> fetches -> raises
if not (size_loc_pass and all(nan_passes)): return ... # checked too late
```

The fetch raises `ValueError: Genomic region out of bounds` inside a
`multiprocessing.Pool` worker, which kills the entire chromosome. Observed on
mm10 chr5 (151,834,684 bp) where a loop caused a fetch to 151,835,000 - 316 bp
past the end - after ~1,364 of 1,423 chunks had already run.

`run_loop_filtering` now returns early when `size_loc_pass` is False, with the
same row shape (`2n+3`) the filters already discard. Loops near chromosome ends
are dropped cleanly instead of aborting the run.

Note there is **no resume**: the intermediate CSV is written every 100 chunks
but never read back, so a restarted chromosome begins from chunk 1.

## Example Workflow

```bash
# 1. Set up environment
mamba activate absloopquantTB

# 2. Calculate P(s) curves
bash 1.1_calculate_P_s_curves_general.sh
cd qshs/250101_calculate_P_s_curves/
sbatch calculate_P_s_curves_pTh17-1.sh
# Wait for jobs to complete...

# 3. Filter loops
bash 1.2_filter_loops.sh
cd qshs/250101_filter_loops_per_chr_fdr0.01/
sbatch filter_loops_pTh17-1_chr1.sh
sbatch filter_loops_pTh17-1_chr2.sh
# ... etc
```