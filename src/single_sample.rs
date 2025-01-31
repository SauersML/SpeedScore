use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

/// Single sample polygenic score calculation.
///
/// `effect_weights` is a map from (chr, pos) -> (effect_allele, effect_weight).
pub fn calculate_polygenic_score(
    path: &str,
    effect_weights: &HashMap<(String, u32), (String, f32)>,
) -> io::Result<(f64, usize, usize, bool)> {
    let file = File::open(path)?;
    let reader = BufReader::with_capacity(1024 * 1024, MultiGzDecoder::new(file)); // 1MB buffer

    // Read entire file lines
    let lines: Vec<String> = reader.lines().collect::<io::Result<_>>()?;

    // Detect whether the VCF uses "chr" prefix by scanning first non‐header line
    let vcf_chr_format = lines.iter()
        .find(|line| !line.starts_with('#'))
        .map(|line| line.starts_with("chr"))
        .unwrap_or(false);

    // We will parallelize over lines, collecting (score, total, matched)
    let (score_sum, total_variants, matched_variants) = lines
        .par_iter()
        .filter(|line| !line.starts_with('#'))
        .map(|line| process_single_sample_line(line, effect_weights))
        .reduce(
            || (0.0, 0, 0),
            |acc, val| (acc.0 + val.0, acc.1 + val.1, acc.2 + val.2),
        );

    Ok((score_sum, total_variants, matched_variants, vcf_chr_format))
}

/// Process a single VCF line for the single‐sample case:
///  - Parse CHR, POS, REF, ALT, sample genotype
///  - If (CHR, POS) in effect_weights, check effect allele vs. REF/ALT
///  - Parse genotype to count effect alleles
/// Returns `(score, total_variants, matched_variants)`.
fn process_single_sample_line(
    line: &str,
    effect_weights: &HashMap<(String, u32), (String, f32)>,
) -> (f64, usize, usize) {
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() < 10 {
        return (0.0, 0, 0); // Malformed line or no genotype
    }

    let chr_raw = parts[0];
    let pos_raw = parts[1];
    let ref_allele = parts[3];
    let alt_allele = parts[4];
    let gt_field = parts[9]; // The sample genotype field (e.g., "0/1", "1/1", "0|1:...")

    // Convert pos to u32
    let pos = match pos_raw.parse::<u32>() {
        Ok(p) => p,
        Err(_) => return (0.0, 0, 0),
    };

    // Normalize chromosome (remove "chr" if present)
    let normalized_chr = chr_raw.trim_start_matches("chr").to_string();

    // If not in effect_weights, skip
    let (effect_allele, weight) = match effect_weights.get(&(normalized_chr.clone(), pos)) {
        Some(x) => x,
        None => return (0.0, 1, 0), // total=1, matched=0
    };

    // Decide if effect_allele is the REF or the ALT. If neither, skip
    let effect_is_ref = effect_allele == ref_allele;
    let effect_is_alt = effect_allele == alt_allele;
    if !effect_is_ref && !effect_is_alt {
        // The scoring file says effect_allele is something else (e.g. "T") 
        // but the VCF has REF="A", ALT="G". No match => skip
        return (0.0, 1, 0);
    }

    // Extract just the genotype itself (e.g. "0/1") from "0/1:..."
    let genotype = gt_field.split(':').next().unwrap_or(".");

    // Count how many effect alleles
    match parse_allele_count(genotype, effect_is_alt) {
        Some(allele_count) => {
            let line_score = *weight as f64 * allele_count as f64;
            (line_score, 1, 1)
        }
        None => {
            // Missing or invalid genotype => skip
            (0.0, 1, 0)
        }
    }
}


/// Helper that counts how many effect alleles are present in `genotype`.
/// If `effect_is_alt` = true, we count `'1'` as effect alleles.
/// If `effect_is_alt` = false, we count `'0'` as effect alleles.
/// Returns None if we see multi‐allelic (e.g. '2') or missing ('.').
fn parse_allele_count(genotype: &str, effect_is_alt: bool) -> Option<u8> {
    let mut count = 0u8;
    for c in genotype.chars() {
        match c {
            '0' if !effect_is_alt => count += 1,
            '1' if effect_is_alt => count += 1,
            '.' | '2' | '3' => return None, // skip multi‐allelic or missing
            '|' | '/' => {} // just a delimiter
            _ => {}
        }
    }
    Some(count)
}
