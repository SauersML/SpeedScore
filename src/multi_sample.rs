use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{self, Read, BufRead, BufReader, Write};
use std::time::Instant;
use std::path::Path;
use flate2::read::MultiGzDecoder;
use indicatif::{ProgressBar, ProgressStyle};

#[derive(Debug)]
pub enum VcfError {
    Io(io::Error),
    Utf8Error(std::string::FromUtf8Error),
}

impl std::fmt::Display for VcfError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            VcfError::Io(err) => write!(f, "I/O error: {}", err),
            VcfError::InvalidFormat(msg) => write!(f, "Invalid format: {}", msg),
            VcfError::Utf8Error(err) => write!(f, "UTF-8 error: {}", err),
        }
    }
}

impl std::error::Error for VcfError {}

impl From<io::Error> for VcfError {
    fn from(error: io::Error) -> Self {
        VcfError::Io(error)
    }
}

impl From<std::string::FromUtf8Error> for VcfError {
    fn from(error: std::string::FromUtf8Error) -> Self {
        VcfError::Utf8Error(error)
    }
}

#[derive(Clone, Default)]
struct SampleData {
    score: f64,
    matched_variants: usize,
    total_variants: usize,
}

fn open_vcf_reader(path: &str) -> Result<BufReader<MultiGzDecoder<File>>, VcfError> {
    let file = File::open(path).map_err(VcfError::Io)?;
    let decoder = MultiGzDecoder::new(file);
    Ok(BufReader::with_capacity(1024 * 1024, decoder)) // 1MB buffer
}

pub fn calculate_polygenic_score_multi(
    vcf_path: &str,
    effect_weights: &HashMap<(String, u32), (String, f32)>,
    output_path: &str,
    debug: bool
) -> Result<(f64, usize, usize, bool), VcfError> {
    let start_time = Instant::now();

    println!("Opening file: {}", vcf_path);
    println!("Effect weights loaded: {} variants", effect_weights.len());

    let mut reader = open_vcf_reader(vcf_path)?;
    let mut header_line = String::new();
    let mut sample_names: Vec<String>;

    // Find the header
    loop {
        reader.read_line(&mut header_line)?;
        if header_line.starts_with("#CHROM") {
            sample_names = header_line.split_whitespace().skip(9).map(String::from).collect();
            break;
        }
        header_line.clear();
    }

    println!("VCF data start found.");
    println!("Sample count: {}", sample_names.len());
    println!("Processing variants...");

    let pb = ProgressBar::new_spinner();
    pb.set_style(ProgressStyle::default_spinner()
        .template("{spinner:.green} [{elapsed_precise}] {msg}")
        .unwrap());
    pb.set_message("Processing...");

    let mut buffer = Vec::new();
    let mut sample_data: Vec<SampleData> = vec![SampleData::default(); sample_names.len()];
    let mut lines_processed = 0;
    let mut last_chr = String::new();
    let mut last_pos = 0;
    let mut vcf_chr_format = false;

    loop {
        buffer.clear();
        let num_lines = reader.read_until(b'\n', &mut buffer)?;
        if num_lines == 0 {
            break;
        }

        lines_processed += 1;

        if !buffer.starts_with(&[b'#']) {
            let result = process_chunk(&buffer, effect_weights, &mut sample_data, debug);
            if let Some((chr, pos, chr_format)) = result {
                if debug && (chr != last_chr || pos > last_pos + 20_000_000) {
                    pb.suspend(|| {
                        println!("\rProcessed up to Chr {}, Pos {:.2}M", chr, pos as f64 / 1_000_000.0);
                        io::stdout().flush().unwrap();
                    });
                    last_chr = chr;
                    last_pos = pos;
                }
                vcf_chr_format = chr_format;
            }
        }

        if lines_processed % 100_000 == 0 {
            let lines_in_k = lines_processed / 1000;
            let variants = sample_data.iter().map(|sd| sd.total_variants).sum::<usize>();
            let matched = sample_data.iter().map(|sd| sd.matched_variants).sum::<usize>();
            pb.set_message(format!(
                "{}K lines, {}K variants, {}K matched",
                lines_in_k,
                variants / 1000,
                matched / 1000
            ));
        }
    }

    pb.finish_with_message("Processing complete");

    let duration = start_time.elapsed();

    write_csv_output(output_path, vcf_path, &sample_names, &sample_data, duration)?;

    let avg_score = sample_data.iter().map(|sd| sd.score).sum::<f64>() / sample_data.len() as f64;
    let total_variants = sample_data.iter().map(|sd| sd.total_variants).sum();
    let matched_variants = sample_data.iter().map(|sd| sd.matched_variants).sum();

    println!("\nFinished processing.");
    println!("Total lines processed: {:.3}K", lines_processed as f64 / 1000.0);
    println!("Results written to: {}", output_path);
    println!("Processing time: {:?}", duration);

    Ok((avg_score, total_variants, matched_variants, vcf_chr_format))
}

/// Processes one chunk of lines (already read from the file).
/// For each line, parse CHR, POS, REF, ALT, then genotypes for each sample.
/// We skip multi‐allelic sites or missing genotypes. 
/// Returns `(last_chr, last_pos, vcf_uses_chr_prefix)`.
fn process_chunk(
    chunk: &[u8],
    effect_weights: &HashMap<(String, u32), (String, f32)>,
    sample_data: &mut [SampleData],
    _debug: bool
) -> Option<(String, u32, bool)> {
    let mut last_chr = String::new();
    let mut last_pos = 0;
    let mut vcf_chr_format = false;

    // Split chunk by newlines
    for line in chunk.split(|&b| b == b'\n') {
        if line.is_empty() || line.starts_with(&[b'#']) {
            continue;
        }

        // Convert line to string
        let line_str = match std::str::from_utf8(line) {
            Ok(s) => s,
            Err(_) => continue, // skip invalid UTF-8
        };

        let parts: Vec<&str> = line_str.split('\t').collect();
        if parts.len() < 10 {
            continue; // skip malformed line
        }

        let chr_raw = parts[0];
        let pos_raw = parts[1];
        let ref_allele = parts[3];
        let alt_allele = parts[4];

        // The 8th column is `FORMAT`; sample genotypes start at index 9
        let genotype_fields = &parts[9..];

        let pos = match pos_raw.parse::<u32>() {
            Ok(p) => p,
            Err(_) => continue,
        };

        last_chr = chr_raw.to_string();
        last_pos = pos;
        vcf_chr_format = chr_raw.starts_with("chr");

        // Normalize chromosome to match how we stored it in effect_weights
        let normalized_chr = chr_raw.trim_start_matches("chr").to_string();

        // If not found in effect_weights, skip
        let (effect_allele, weight) = match effect_weights.get(&(normalized_chr.clone(), pos)) {
            Some(x) => x,
            None => {
                // Still count total_variants for each sample?
                for sample in sample_data.iter_mut() {
                    sample.total_variants += 1;
                }
                continue;
            }
        };

        // Check if effect allele is REF or ALT. Otherwise skip
        let effect_is_ref = effect_allele == ref_allele;
        let effect_is_alt = effect_allele == alt_allele;
        if !effect_is_ref && !effect_is_alt {
            // Increase total_variants but not matched
            for sample in sample_data.iter_mut() {
                sample.total_variants += 1;
            }
            continue;
        }

        // At this point, we have a matched variant that matters for scoring
        // Increase total_variants and matched_variants for each sample
        for (sample, genotype_field) in sample_data.iter_mut().zip(genotype_fields) {
            sample.total_variants += 1;
            sample.matched_variants += 1;

            // The genotype might look like "0/1:..." so we isolate the GT
            let gt = genotype_field.split(':').next().unwrap_or(".");

            // Count how many effect alleles
            match parse_allele_count(gt, effect_is_alt) {
                Some(allele_count) => {
                    sample.score += (*weight as f64) * (allele_count as f64);
                }
                None => {
                    // skip if missing or multi-allelic
                }
            }
        }
    }

    Some((last_chr, last_pos, vcf_chr_format))
}

/// Identical to the single-sample helper (move to common later):
/// If `effect_is_alt`, we count '1' as effect alleles; otherwise we count '0'.
fn parse_allele_count(gt: &str, effect_is_alt: bool) -> Option<u8> {
    let mut count = 0u8;
    for c in gt.chars() {
        match c {
            '0' if !effect_is_alt => count += 1,
            '1' if effect_is_alt => count += 1,
            '.' | '2' | '3' => return None, // skip multi-allelic or missing
            '|' | '/' => {}
            _ => {}
        }
    }
    Some(count)
}

fn write_csv_output(
    output_path: &str,
    vcf_path: &str,
    sample_names: &[String],
    sample_data: &[SampleData],
    duration: std::time::Duration
) -> Result<(), VcfError> {
    let path = Path::new(output_path);
    let prefix = path.parent().unwrap_or_else(|| Path::new(""));
    std::fs::create_dir_all(prefix).map_err(VcfError::Io)?;

    let mut file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(output_path)
        .map_err(VcfError::Io)?;

    writeln!(file, "VCF_File,Sample_Name,Polygenic_Score,Calculation_Time_Seconds,Total_Variants,Matched_Variants")
        .map_err(VcfError::Io)?;

    for (name, data) in sample_names.iter().zip(sample_data.iter()) {
        writeln!(
            file,
            "{},{},{:.6},{:.6},{},{}",
            vcf_path,
            name,
            data.score,
            duration.as_secs_f64(),
            data.total_variants,
            data.matched_variants
        ).map_err(VcfError::Io)?;
    }

    Ok(())
}
