# Hash Calculator

A fast, parallel hash calculator for single-cell RNA-seq data stored in h5ad format.

## Description

This tool computes SHA256 hashes for individual cells in single-cell datasets based on a selected subset of genes. It reads AnnData objects from h5ad files, filters the gene expression matrix to a specified gene list, and generates reproducible hashes for each cell.

## Installation

Build from source using Cargo:

```bash
cargo build --release
```

The binary will be located at `target/release/hash_calculator`.

## Usage

```bash
hash_calculator --input <INPUT_FILE> --intersect <GENE_LIST> --output <OUTPUT_FILE> --jobs <NUM_THREADS> --seed <SEED>
```

### Arguments

- `-i, --input <INPUT_FILE>` - Path to input h5ad file containing single-cell data
- `--intersect <GENE_LIST>` - Path to text file containing gene names (one per line) to include in hash calculation
- `-o, --output <OUTPUT_FILE>` - Path to output CSV file for cell hashes
- `-j, --jobs <NUM_THREADS>` - Number of parallel threads to use
- `-s, --seed <SEED>` - Random seed for reproducibility (integer)

### Example

```bash
hash_calculator \
  --input data.h5ad \
  --intersect genes.txt \
  --output cell_hashes.csv \
  --jobs 8 \
  --seed 42
```

## Output Format

The tool generates a CSV file with two columns:

```
cell_name,hash
cell_001,a3f5b1c2d4e6f7a8b9c0d1e2f3a4b5c6d7e8f9a0b1c2d3e4f5a6b7c8d9e0f1a2
cell_002,b4f6c2d3e5f7a8b9c0d1e2f3a4b5c6d7e8f9a0b1c2d3e4f5a6b7c8d9e0f1a2b3
...
```

## Features

- Supports sparse CSR matrices with multiple numeric types (i8, i16, i32, i64, u8, u16, u32, u64, f32, f64)
- Parallel processing using Rayon for improved performance
- Memory-efficient chunked processing of large datasets
- SHA256 hashing for cryptographic-quality fingerprints

## Requirements

- Rust 1.70 or later
- Input data must be in h5ad format (AnnData HDF5)
