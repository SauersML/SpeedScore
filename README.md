# SpeedScore

## Overview

This Rust-based tool calculates polygenic scores from VCF (Variant Call Format) files using a provided scoring file.

## Features

- Fast polygenic score calculation
- Parallel processing
- Memory-mapped file I/O for efficient data handling
- Supports large VCF files and scoring files

## Prerequisites

- Rust (latest stable version)
- Cargo (Rust's package manager)

## Installation

1. Clone the repository:
   ```
   git clone https://github.com/ScottSauers/SpeedScore.git
   cd SpeedScore
   ```

2. Build the project:
   ```
   cargo build --release
   ```

## Usage

Run the program with the following command:

```
cargo run --release -- -v <path_to_vcf_file> -s <path_to_scoring_file> --output <path_to_output_file> [--info]
```

### Command-line Arguments

- `-v, --vcf <FILE>`: Path to the input VCF file
- `-s, --scoring <FILE>`: Path to the scoring file
- `--output <FILE>`: Path to the output file
- `--info`: (Optional) Display detailed information about the calculation

### Example

```
cargo run --release -- -v /path/to/your/file.vcf -s /path/to/your/scoring.txt --output /path/to/your/output.txt --info
```

## File Formats

### VCF File
The input should be a standard VCF file. The tool expects the chromosome, position, and genotype information.

### Scoring File
The scoring file is expected to be in PGS Catalog format; that is, a tab-separated file with the following columns:
1. Chromosome
2. Position
3. Effect allele
4. Other allele
5. Effect weight

Example:
```
1   760912  C   T   8.06914e-05
1   846808  C   T   0.000365935
1   861808  A   G   -0.000241058
```

## Output

When the `--info` flag is used, additional information is displayed in the console. The tool can generate a tab-separated output file containing:
- VCF file path
- Scoring file path
- Calculated polygenic score
- Calculation time
- Total variants processed
- Matched variants
- Number of variants in the scoring file
