# LD Adjuster

A high-performance C++ implementation for linkage disequilibrium (LD) adjustment of genomic relationship matrices, optimized for large-scale genetic data.

## Features

- **High Performance**: Optimized for dense matrices up to 50k×50k
- **Memory Efficient**: Uses single-precision floating point arithmetic
- **Robust**: Handles missing data and invariant SNPs automatically
- **Compatible**: Reads standard PLINK raw files and GRM binary format
- **Fast**: Uses Eigen library with optional BLAS acceleration

## Requirements

### Dependencies
- **Eigen3**: Linear algebra library (≥3.3)
- **C++ Compiler**: Supporting C++11 (GCC ≥4.8, Clang ≥3.3)
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

# Debug build (for development)
make debug

# Clean build files
make clean

# Show build configuration
make info
```

## Usage

### Command Line Interface
```bash
./ld_adjuster --raw <input.raw> --grm <grm_prefix> --h2 <heritability> --out <output_prefix> [--threads N]
```

### Parameters
- `--raw`: PLINK raw format genotype file
- `--grm`: Prefix for GRM files (expects .grm.iid and .grm.bin)
- `--h2`: Heritability value (0.0 to 1.0)
- `--out`: Output file prefix
- `--threads N`: Number of threads to use (optional, default: all available)
- `--num_threads N`: Alias for --threads
- `--help`: Show usage information

### Example
```bash
./ld_adjuster --raw genotypes.raw --grm population_grm --h2 0.5 --out results --threads 8
```

## Input File Formats

### PLINK Raw File (.raw)
Tab-delimited text file with:
- Header row: `FID IID PAT MAT SEX PHENOTYPE SNP1 SNP2 ...`
- Data rows: Individual information + genotypes (0/1/2)

### GRM Files
- **`.grm.iid`**: Text file with individual IDs (FID IID per line)
- **`.grm.bin`**: Binary file containing:
  - Number of individuals (int32)
  - Number of markers (float32)
  - Lower triangular matrix values (float32)

## Output Files

The program generates three output files:

### 1. Summary File (`.neff`)
```
Effective sample size: 1234
Number of matched individuals: 2000
Number of SNPs: 45000
Heritability: 0.50
Computation time (ms): 15432
Peak memory usage (MB): 8240
```

### 2. LD Matrix (`.ldmatrix`)
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
rs1234567
rs2345678
rs3456789
...
```

## Algorithm Details

The LD adjustment algorithm follows these steps:

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
For a matrix with `n` SNPs and `m` individuals:
- **GRM matrix**: m × m × 4 bytes
- **SNP matrix**: n × n × 4 bytes  
- **Raw data**: m × n × 4 bytes
- **Total**: ~(m² + n² + mn) × 4 bytes

**Examples:**
- **50k SNPs, 10k individuals**: ~12 GB RAM
- **25k SNPs, 5k individuals**: ~2.7 GB RAM
- **10k SNPs, 2k individuals**: ~0.4 GB RAM

### Computation Time
With optimized BLAS (approximate times):
- **10k SNPs**: 1-5 minutes
- **25k SNPs**: 15-45 minutes  
- **50k SNPs**: 1-3 hours

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

### Getting Help

For questions or issues:
1. Check this README for common solutions
2. Verify input file formats match specifications
3. Test with smaller datasets first
4. Monitor memory usage during execution

## Technical Details

- **Precision**: Single-precision floating point (32-bit)
- **Matrix Storage**: Dense format in memory
- **Decomposition**: LLT (Cholesky) for positive definite matrices
- **Missing Data**: Replaced with SNP mean before standardization
- **Threading**: Depends on underlying BLAS implementation

## License

This implementation is provided for academic and research use.