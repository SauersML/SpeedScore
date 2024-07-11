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
cargo run --release -- -v <path_to_vcf_file> -s <path_to_scoring_file> --output <path_to_output_file>
```
Add the ```--info``` flag for additional information.

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

## Multi-sample VCF
SpeedScore also supports multi-sample VCFs. For example, it can run a VCF containing the 1000 Genomes dataset.

A multi-sample VCF will have the following output format:

| VCF_File                    | Sample_Name | Polygenic_Score | Calculation_Time_Seconds | Total_Variants | Matched_Variants |
|-----------------------------|-------------|-----------------|--------------------------|----------------|------------------|
| /path/to/your/file.vcf.gz   | HG00096     | -0.002532       | 1164.698128              | 45693          | 45693            |
| /path/to/your/file.vcf.gz   | HG00097     | 0.007797        | 1164.698128              | 45693          | 45693            |
| /path/to/your/file.vcf.gz   | HG00099     | -0.009080       | 1164.698128              | 45693          | 45693            |
| /path/to/your/file.vcf.gz   | HG00100     | -0.012225       | 1164.698128              | 45693          | 45693            |



Here is an example of how this might be done:

```
# Clone the repository
git clone https://github.com/ScottSauers/SpeedScore.git
cd SpeedScore
sudo apt-get update
sudo apt-get install aria2

# Download your file
aria2c -x 16 -s 16 -o 1000G_file.tar.zst https://your/url/to/1000genomes/1000G_file.tar.zst
zstd -d 1000G_file.tar.zst
tar -xf 1000G_file.tar
sudo apt install plink2

# We would like the file in VCF format. Find the relevant file from the decompressed binary
plink2 --pfile GRCh38_1000G_ALL 'vzs' --allow-extra-chr --chr 1-22,X,Y,MT \
  --export vcf bgz \
  --out GRCh38_1000G_ALL_vcf

# Download a scoring file of your choice. For example, PGS003725
curl -O https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS003725/ScoringFiles/Harmonized/PGS003725_hmPOS_GRCh38.txt.gz
gunzip PGS003725_hmPOS_GRCh38.txt.gz
sudo apt  install cargo -y
cargo build --release

# Run the command with --info flag for extra information
cargo run --release -- -v /home/your/path/to/SpeedScore/GRCh38_1000G_ALL_vcf.vcf.gz -s /home/your/path/to/SpeedScore/PGS003725_hmPOS_GRCh38.txt --output /home/your/path/to/SpeedScore/output_multi_sample.txt --info
```
For a scoring with 1,296,172 variants, it took 1165.299863 seconds, or (3202/1165.299863) = ~2.75 seconds per ~146 million variant file, on a machine with these specs:

- **OS**: Linux, Ubuntu 22.04.1
- **Architecture**: x86_64
- **Processor**: Intel(R) Xeon(R) Platinum 8370C CPU @ 2.80GHz
- **Cores**: 4 (2 cores per socket, 2 threads per core)
- **CPU max MHz**: 2800.0000
- **CPU min MHz**: 800.0000
- **L1 Cache**: 96 KiB (2 instances)
- **L2 Cache**: 2.5 MiB (2 instances)
- **L3 Cache**: 48 MiB (1 instance)
- **Available Memory**: 30 GiB
- **Swap**: 0B
- **Storage**: 68GB NVMe disk
