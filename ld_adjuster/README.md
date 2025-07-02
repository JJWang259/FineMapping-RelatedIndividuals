# LD Adjuster

A C++ program for adjusting linkage disequilibrium (LD) matrix with genomic relationship matrix.

## Features

- Handles missing data and invariant SNPs automatically
- Reads PLINK raw files and [MPH](https://jiang18.github.io/mph/) GRM files
- Uses Eigen library with optional BLAS acceleration

## Requirements

### Dependencies
- **Eigen3**: Linear algebra library (≥3.3)
- **C++ Compiler**: Supporting C++11 (GCC ≥4.8 or Clang ≥3.3)
- **Optional**: OpenBLAS or Intel MKL for maximum performance

### Installation on Linux (Ubuntu/Debian)
```bash
# Install dependencies
sudo apt-get update
sudo apt-get install build-essential libeigen3-dev

# For optimal performance, install OpenBLAS
sudo apt-get install libopenblas-dev

# Clone or download the source files
# ld_adjuster.cpp and Makefile
```

### Installation on Red Hat/CentOS
```bash
# Install dependencies
sudo yum install gcc-c++ eigen3-devel

# For RHEL 8+/CentOS 8+
sudo dnf install gcc-c++ eigen3-devel openblas-devel
```

## Building

```bash
# Standard build
make

# Clean build files
make clean

# Show build configuration
make info
```

## Pre-compiled Binary
https://github.com/JJWang259/FineMapping-RelatedIndividuals/releases/latest

## Usage

### Command Line Interface
```bash
./ld_adjuster --raw <input.raw> --grm <grm_prefix> --h2 <heritability> --out <output_prefix> [--threads N]
```

### Parameters
- `--raw`: PLINK raw format genotype file
- `--grm`: Prefix for [MPH](https://jiang18.github.io/mph/) GRM files (expects .grm.iid and .grm.bin)
- `--h2`: Heritability value (0.0 to 1.0)
- `--out`: Output file prefix
- `--threads N`: Number of threads to use (optional, default: all available)
- `--num_threads N`: Alias for --threads
- `--help`: Show usage information

### Example
```bash
./ld_adjuster --raw test.raw --grm test --h2 0.5 --out test --threads 8
```

## Input File Formats

### PLINK Raw File (.raw)
Space-delimited text file with:
- Header row: `FID IID PAT MAT SEX PHENOTYPE SNP1 SNP2 ...`
- Data rows: Individual information + genotypes (0/1/2)

### [MPH](https://jiang18.github.io/mph/) GRM Files
- **`.grm.iid`**: Text file with individual IDs (IID per line)
- **`.grm.bin`**: Binary file containing:
  - Number of individuals (int32)
  - Number of markers or sum of weights (float32)
  - Lower triangular matrix values (float32)

## Output Files

The program generates three output files:

### 1. Summary File (`.summary`)
```
Effective sample size: 1234
Number of matched individuals: 2000
Number of SNPs: 5000
Heritability: 0.50
Computation time (ms): 432
Peak memory usage (MB): 41
```

### 2. LD Matrix (`.ld`)
Space-delimited correlation matrix (SNPs × SNPs):
```
1.000000 0.234567 -0.123456 ...
0.234567 1.000000 0.345678 ...
-0.123456 0.345678 1.000000 ...
...
```

### 3. SNP IDs (`.snpids`)
One SNP identifier per line:
```
rs1234567_T
rs2345678_A
rs3456789_G
...
```

## Algorithm Details

The LD matrix adjustment algorithm follows these steps:

1. **Data Matching**: Finds intersection of individuals between raw and GRM files
2. **Quality Control**: Removes SNPs with zero variance
3. **Heritability Adjustment**: Modifies GRM diagonal: `G_ii = η × G_ii + 1` where `η = h²/(1-h²)`
4. **Standardization**: Centers and scales genotype data
5. **Cholesky Decomposition**: Computes `L` such that `G = LL^T`
6. **Linear System**: Solves `LX* = X` for adjusted genotypes
7. **Correlation Matrix**: Computes `R* = (X*)^T X*`
8. **Normalization**: Standardizes to correlation matrix

## Performance Optimization

### Memory Requirements
For a matrix with `m` SNPs and `n` individuals:
- **GRM matrix**: n × n × 4 bytes
- **SNP matrix**: m × m × 4 bytes  
- **Raw data**: n × m × 4 bytes
- **Total**: ~(m² + 2n² + mn) × 4 bytes

### Optimization Tips
1. **Use OpenBLAS**: Install `libopenblas-dev` for best performance
2. **Sufficient RAM**: Ensure adequate memory to avoid swapping
3. **Fast Storage**: Use SSD for large input files
4. **Compiler Flags**: The Makefile uses `-O3 -march=native` for optimization

## Troubleshooting

### Common Issues

**"Eigen3 not found"**
```bash
# Install Eigen3 development headers
sudo apt-get install libeigen3-dev
```

**"Cholesky decomposition failed"**
- Check GRM matrix is positive definite
- Verify heritability parameter is reasonable (0.1-0.9)
- Ensure sufficient matched individuals

**"Out of memory"**
- Reduce problem size or use machine with more RAM
- Check if swap space is configured

**Slow performance**
- Install optimized BLAS: `sudo apt-get install libopenblas-dev`
- Rebuild after installing BLAS: `make clean && make`
- Use multiple threads such as `--threads 8`

### Getting Help

For questions or issues:
1. Check this README for common solutions
2. Verify input file formats match specifications
3. Monitor memory usage during execution

## Technical Details

- **Precision**: Single-precision floating point (32-bit)
- **Matrix Storage**: Dense format in memory
- **Decomposition**: LLT (Cholesky) for positive definite matrices
- **Missing Data**: Replaced with SNP mean before standardization
- **Threading**: Depends on underlying BLAS implementation

## License

This implementation is provided for academic and research use.
