# AbsLoopQuant_TB - Absolute Loop Quantification Toolkit

This package provides tools for calculating P(s) curves from Hi-C data and filtering loops based on quantitative criteria. It implements the AbsLoopQuant approach for identifying high-confidence chromatin loops.

## Overview

The AbsLoopQuant workflow consists of four steps:

1. **Calculate P(s) curves** - Compute probability of contact as a function of genomic distance from Hi-C cooler files
2. **Filter loops** - Filter loop calls based on size, location, NaN regions, read counts, and global maximum distance criteria
3. **Concatenate** - Join the per-chromosome output of step 2 into one loop table per condition
4. **Intersect** - Ask which of those loops are also called by fithic in each individual replicate

## Files

- **`1.1_calculate_P_s_curves_general.py`** - Calculate P(s) curves from Hi-C cooler files
- **`1.1_calculate_P_s_curves_general.sh`** - Bash script to generate SLURM jobs for P(s) curve calculation
- **`1.2_filter_loops.py`** - Filter loops per chromosome using quantitative criteria
- **`1.2_filter_loops.sh`** - Bash script to generate SLURM jobs for loop filtering
- **`1.3_concat_loops.sh`** - Concatenate the per-chromosome loop tables into one per condition
- **`1.4_intersect_replicates.py`** - Intersect a loop set with per-replicate fithic calls
- **`1.4_intersect_replicates.sh`** - Driver for the intersection across conditions and replicates
- **`looptools.py`** - Helper module with loop analysis utilities
- **`0.0_create_absloopquantTB_env.sh`** - Builds the `absloopquantTB` mamba environment and verifies its imports
- **`absloopquantTB_env.yml`** - Conda/mamba environment spec

## Setup

Requires `mamba` (or `conda`). Everything installs from `conda-forge` and
`bioconda`.

### 1. Create the environment

```bash
bash 0.0_create_absloopquantTB_env.sh
```

This solves `absloopquantTB_env.yml`, then verifies that `numpy`, `pandas`,
`scipy`, `cooler`, `cooltools`, `cv2`, `matplotlib` and `looptools` all import -
a successful solve alone is not proof the environment works.

Two pins in the spec are not cosmetic, and relaxing either one breaks the
pipeline at import rather than at runtime:

- **`numpy<2`** - `cooler`, `cooltools`, `py-opencv` and `numba` ship C
  extensions built against the numpy 1.x ABI, and numpy 2 breaks them at import.
- **`matplotlib-base<3.9`** - `cooltools.lib.plotting` imports `register_cmap`,
  removed in matplotlib 3.9, and its fallback targets a submodule that has never
  existed. `looptools` imports that module at load time, so **step 1 fails at
  import** without this pin.

Options:

```bash
ENV_NAME=absloopquantTB2 bash 0.0_create_absloopquantTB_env.sh   # build side-by-side
ENV_PREFIX=/path/to/envs/absloopquantTB bash 0.0_create_absloopquantTB_env.sh
```

`ENV_PREFIX` builds the environment at a path you choose instead of inside the
conda root - useful when your home directory is small or quota'd.

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

### Step 3: Concatenate Loops Across Chromosomes

Step 2 writes one file per chromosome. Step 3 joins them into a single table per
condition, and renames the columns to fithic's, which is what lets step 4 join
the two on coordinates.

```bash
bash 1.3_concat_loops.sh
```

Unlike steps 1 and 2 this does the work directly rather than generating SLURM
scripts - it is a `cat`/`awk` pass over ~20 small files.

| | |
|---|---|
| Input | `<results>/filtered_loops_per_chr_fdr<fdr>/{condition}/filtered_loops.fdr<fdr>.<chrom>.txt` |
| Output | `<results>/filtered_loops_<res>kb_fdr<fdr>/{condition}/{condition}.coords.fdr<fdr>.txt` |
| Columns | `chr1  fragmentMid1  chr2  fragmentMid2  size` |

The `chroms` array near the top lists mouse autosomes (`chr1`..`chr19`). Change
it for other assemblies - human needs `chr1`..`chr22` - and it must match what
step 2 was actually run on.

An existing output file is left alone; pass `FORCE=1` to rebuild it. A missing
chromosome is reported and skipped rather than aborting the condition, so check
the `N chromosomes, M missing` line before using the result.

### Step 4: Intersect With Per-Replicate fithic Calls

Steps 1-3 produce high-confidence loops from the combined sample. Step 4 asks,
for each loop, which individual replicates fithic also called it in.

```bash
bash 1.4_intersect_replicates.sh
```

The fithic files are chosen in this order of precedence:

1. `FITHIC_FILES` - a space-separated list, applied to every condition:
   ```bash
   FITHIC_FILES="/path/rep1.significances.txt /path/rep2.significances.txt.gz" \
     bash 1.4_intersect_replicates.sh
   ```
2. the `fithicFiles` array in the script, if you fill it in
3. otherwise derived from `FITHIC_DIR` and `replicateNamesList` using the layout
   `$FITHIC_DIR/<replicate>/fithic/<resolution>/<replicate>.<template>`

`.gz` files are read directly - there is no decompression step.

| | |
|---|---|
| Input | step 3's `{condition}.coords.fdr<fdr>.txt`, plus one fithic file per replicate |
| Output | `<results>/fithic_filtered_loops_bioreplicates_<res>kb_fdr<fdr>/loops/{condition}/<replicate>.<...>.fdr<fdr>.significances.txt` |

Each output is an inner join on `chr1`, `fragmentMid1`, `chr2`, `fragmentMid2`,
carrying the replicate's own `contactCount`, `p-value`, `q-value` and biases
through. The run prints how many of the loop set each replicate recovered.

The replicate table is **not** q-value filtered by default: the loop set is
already FDR-filtered, and the question here is presence, not independent
significance. Pass `--apply_fdr` to the Python script if you want both.

To run the intersection on its own, outside the driver:

```bash
python3 1.4_intersect_replicates.py \
    --combined_loops_file /path/to/{condition}.coords.fdr0.01.txt \
    --replicate_files rep1.significances.txt rep2.significances.txt.gz \
    --output_dir /path/to/output \
    --fdr_threshold 0.01 \
    --verbose
```

## Configuring the job generators

`1.1_calculate_P_s_curves_general.sh` and `1.2_filter_loops.sh` resolve this
repo from their own location, so a clone works from any path. The data project
is site-specific and is set by environment variable rather than by editing the
files:

```bash
BASE_DIR=/path/to/your/hic-project bash 1.2_filter_loops.sh
WORKING_DIR=/path/to/AbsLoopQuant_TB bash 1.2_filter_loops.sh   # code elsewhere
```

| Variable | Used by | Meaning |
|---|---|---|
| `WORKING_DIR` | all | This repo. Defaults to the script's own directory. |
| `BASE_DIR` | all | The data project root. |
| `RESULTS_DIR` | 1.3, 1.4 | Where step 2 wrote its output. Set it to keep steps 2-4 pointed at one place. |
| `FITHIC_DIR` | 1.4 | Root of per-replicate fithic output, when deriving paths. |
| `FITHIC_FILES` | 1.4 | Explicit space-separated list of fithic files. |
| `FORCE` | 1.3 | Overwrite an existing concatenated table. |

`BASE_DIR` is expected to hold HiC-Pro-style output:

```
$BASE_DIR/<per-replicate matrix dir>/{condition}-<rep>/cool/{condition}-<rep>.mcool
$BASE_DIR/<combined matrix dir>/{condition}/cool/{condition}.mcool
$BASE_DIR/<combined matrix dir>/{condition}/fithic/<res>/{condition}.*.significances.txt.gz
```

The `subsets`, `resolution`, `nproc` and `fdrThreshold` values near the top of
each generator, and the `perReplicateDir` / `combinedReplicateDir` layout, need
editing to match your project. The generators require SLURM (`sbatch`); to run
without a scheduler, call the two Python scripts directly as shown above.

## Example Workflow

`{condition}` below is one entry of the `subsets` array in the generators, and
`{date}` is the `YYMMDD` stamp the generators put on their output directory.

```bash
# 1. Set up environment
mamba activate absloopquantTB

# 2. Calculate P(s) curves
bash 1.1_calculate_P_s_curves_general.sh
cd qshs/{date}_calculate_P_s_curves/
sbatch calculate_P_s_curves_{condition}.sh
# Wait for jobs to complete...

# 3. Filter loops
bash 1.2_filter_loops.sh
cd qshs/{date}_filter_loops_per_chr_fdr0.01/
sbatch filter_loops_{condition}_chr1.sh
sbatch filter_loops_{condition}_chr2.sh
# ... etc
# Wait for all chromosomes to complete...

# 4. Concatenate the per-chromosome tables
bash 1.3_concat_loops.sh

# 5. Intersect with the per-replicate fithic calls
bash 1.4_intersect_replicates.sh
```

Each step consumes the previous one's output, so they are strictly ordered.
Two places where that bites:

- Step 2 must not start until step 1 has finished for **every** sample it will
  read - the filter reads the P(s) curve of each replicate and of the combined
  sample.
- Step 3 must not start until every chromosome job from step 2 has finished. It
  will happily concatenate a partial set and only warn about what is missing.