# LD Adjuster by Jicai Jiang

A C++ program for adjusting linkage disequilibrium (LD) matrix with genomic relationship matrix and optional covariates.

## Features

- Handles missing data and invariant SNPs automatically
- Reads PLINK raw files and [MPH](https://jiang18.github.io/mph/) GRM files
- Uses Eigen library with optional MKL acceleration
- Supports individual-specific error weights
- Includes covariate adjustment for confounding factors

## Compilation

### Dependencies
- [**Eigen 3**](https://eigen.tuxfamily.org/index.php?title=Main_Page): Linear algebra library
- [**Intel oneAPI HPC Toolkit**](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html): The toolkit offers maximum performance. 

```bash
icpx -Wall -O3 -qmkl -qopenmp -std=c++11 -DEIGEN_USE_MKL_ALL -I/path/to/eigen-3.4.0 -o ld_adjuster ld_adjuster.cpp
```

## Binary Executable
https://github.com/JJWang259/FineMapping-RelatedIndividuals/releases/latest

## Usage

### Command Line Interface
```bash
./ld_adjuster --raw <input.raw> --grm <grm_prefix> --h2 <heritability> --out <output_prefix> [options]
```

### Parameters
- `--raw`: PLINK raw format genotype file
- `--grm`: Prefix for [MPH](https://jiang18.github.io/mph/) GRM files (expects .grm.iid and .grm.bin)
- `--h2`: Heritability value (0.0 to 1.0)
- `--out`: Output file prefix
- `--threads N`: Number of threads to use (optional, default: all available)
- `--num_threads N`: Alias for --threads
- `--error_weight_file`: CSV file with individual error weights (requires --error_weight_name)
- `--error_weight_name`: Column name for error weights (requires --error_weight_file)
- `--covariate_file`: CSV file with covariates (requires --covariate_names)
- `--covariate_names`: Comma-separated covariate names or 'all' (requires --covariate_file)
- `--help`: Show usage information

### Examples
```bash
# Basic usage
./ld_adjuster --raw test.raw --grm test --h2 0.5 --out test --threads 8

# With individual error weights
./ld_adjuster --raw test.raw --grm test --h2 0.5 --out test \
  --error_weight_file weights.csv --error_weight_name weight_col

# With covariate adjustment
./ld_adjuster --raw test.raw --grm test --h2 0.5 --out test \
  --covariate_file covariates.csv --covariate_names age,sex,pc1,pc2

# With covariate adjustment (all columns)
./ld_adjuster --raw test.raw --grm test --h2 0.5 --out test \
  --covariate_file covariates.csv --covariate_names all

# Complete analysis with all features
./ld_adjuster --raw test.raw --grm test --h2 0.5 --out test \
  --error_weight_file weights.csv --error_weight_name weight_col \
  --covariate_file covariates.csv --covariate_names age,sex,pc1 \
  --threads 8
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

### Error Weight File (.csv) - Optional
CSV file with individual-specific error weights:
- Header row with column names
- First column: Individual ID (IID) matching other files
- Subsequent columns: Named weight values for different traits
- Example:
```
IID,trait1_weight,trait2_weight
SAMPLE001,1.2,0.8
SAMPLE002,0.9,1.1
```

### Covariate File (.csv) - Optional
CSV file with covariates for adjustment:
- Header row with column names
- First column: Individual ID (IID) matching other files
- Subsequent columns: Named covariate values
- Example:
```
IID,age,sex,pc1,pc2
SAMPLE001,45,1,0.12,-0.05
SAMPLE002,32,0,-0.08,0.23
```

## Output Files

The program generates three output files:

### 1. Summary File (`.summary`)
```
Effective sample size: 1234
Number of matched individuals: 2000
Number of SNPs: 5000
Heritability: 0.50
Number of covariates: 4
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

1. **Data Matching**: Finds intersection of individuals across all specified files (genotype, GRM, and optional weight/covariate files)
2. **Quality Control**: Removes SNPs with zero variance
3. **Error Weight Adjustment**: Modifies GRM diagonal: `G_ii = η × G_ii + R_ii` where `η = h²/(1-h²)` and `R_ii` are individual error weights (default: 1.0)
4. **Standardization**: Centers and scales genotype and covariate data
5. **Cholesky Decomposition**: Computes `L` such that `G = LL^T`
6. **Linear System**: Solves `LX* = X` for adjusted genotypes and `LQ* = Q` for adjusted covariates (if present)
7. **Correlation Matrix**: 
   - Without covariates: `R* = (X*)^T X*`
   - With covariates: `R* = (X*)^T X* - (X*)^T Q* ((Q*)^T Q*)^(-1) (Q*)^T X*`
8. **Normalization**: Standardizes to correlation matrix

## Performance Optimization

### Memory Requirements
For a matrix with `m` SNPs, `n` individuals, and `c` covariates:
- **GRM matrix**: n × n × 4 bytes
- **SNP matrix**: m × m × 4 bytes  
- **Raw data**: n × m × 4 bytes
- **Covariate data**: n × c × 4 bytes
- **Total**: ~(m² + 2n² + mn + nc) × 4 bytes

### Optimization Tips
1. **Use Intel oneAPI**: Install [Intel oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html) for best performance
2. **Sufficient RAM**: Ensure adequate memory to avoid swapping
3. **Fast Storage**: Use SSD for large input files
4. **Compiler Flags**: The Makefile uses `-O3 -march=native` for optimization

## Troubleshooting

### Common Issues

**"Eigen3 not found"**
- Verify the path in `-I/path/to/eigen-3.4.0` matches your actual Eigen installation directory.

**"Cholesky decomposition failed"**
- Check GRM matrix is positive definite
- Verify heritability parameter is reasonable
- Ensure sufficient matched individuals

**"Covariate matrix is not full column rank"**
- Remove linearly dependent covariates
- Check for constant columns or perfect multicollinearity
- Ensure sufficient individuals relative to number of covariates

**"Column 'X' not found in CSV file"**
- Verify column names match exactly (case-sensitive)
- Check CSV file has proper header row
- Ensure no extra whitespace in column names

**"--error_weight_file requires --error_weight_name"**
- Both error weight options must be used together
- Specify both the CSV file and the column name to use

**"--covariate_file requires --covariate_names"**
- Both covariate options must be used together
- Specify both the CSV file and the column names (or 'all')

**"Out of memory"**
- Reduce problem size or use machine with more RAM
- Check if swap space is configured

**Slow performance**
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
- **Missing Data**: 
  - **SNPs**: Replaced with SNP mean before standardization
  - **Covariates/Error Weights**: Individuals with any missing values are excluded from analysis

## License

This implementation is provided for academic and research use.
