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
- **`create_absloopquantTB_env.sh`** - Script to export the `absloopquantTB` mamba environment
- **`absloopquantTB_env.yml`** - Conda/mamba environment file with all dependencies

## Setup

### 1. Create the environment from the yml file:

```bash
mamba env create -n absloopquantTB -f absloopquantTB_env.yml
```

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
    --P_s_curves_dir /path/to/P_s_curves \
    --looptools_path /path/to/looptools \
    --resolution 10000 \
    --verbose
```

#### Arguments

**Required:**
- `--per_replicate_files`: List of per-replicate `.mcool` file paths
- `--combined_replicate_file`: Path to combined replicate `.mcool` file
- `--per_replicate_names`: List of per-replicate names
- `--combined_replicate_name`: Name of combined replicate
- `--P_s_curves_dir`: Output directory for P(s) curves

**Optional:**
- `--looptools_path`: Path to looptools directory (default: `/home/bbabatunde/packages/25-09-absloopquant/AbsLoopQuant_analysis_code`)
- `--resolution`: Resolution in base pairs (default: 10000)
- `--verbose`: Enable verbose output (default: True)

#### Output

P(s) curve files are saved as:
```
{P_s_curves_dir}/{sample_name}/P_s_{resolution}bp.txt
```

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