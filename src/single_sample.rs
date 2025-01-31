use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::cell::RefCell;
use std::rc::Rc;
use crate::common::ChromosomeFormat;

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

fn detect_vcf_chr_format(lines: &[String]) -> bool {
    lines.iter()
        .find(|line| !line.starts_with('#'))
        .map(|line| line.split_once('\t').unwrap().0.starts_with("chr"))
        .unwrap_or(false)
}

fn process_line(
    line: &str,
    effect_weights: &HashMap<(String, u32), f32>,
    index: usize,
) -> (f64, usize, usize) {
    if index < 5 {
        println!("Raw VCF line {}: {}", index + 1, line);
    }

    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() < 10 {
        if index < 5 {
            println!("Invalid VCF line format: {:?}", parts);
        }
        return (0.0, 1, 0);
    }

    let chr = parts[0];
    let normalized_chr = chr.trim_start_matches("chr").to_string();

    if let Ok(pos) = parts[1].parse::<u32>() {
        if index < 5 {
            println!("Processing variant (example): chr={}, pos={:?}", chr, pos);
        }

        if let Some(&weight) = effect_weights.get(&(normalized_chr, pos)) {
            let genotype = parts[9];
            let allele_count = match genotype.chars().next() {
                Some('0') => 0,
                Some('1') => {
                    if genotype.chars().nth(2) == Some('1') {
                        2
                    } else {
                        1
                    }
                }
                _ => return (0.0, 1, 0),
            };

            let score = f64::from(weight) * allele_count as f64;

            if index < 5 {
                println!(
                    "Matched variant: chr={}, pos={:?}, weight={}, allele_count={}",
                    chr, pos, weight, allele_count
                );
            }

            (score, 1, 1)
        } else {
            (0.0, 1, 0)
        }
    } else {
        if index < 5 {
            println!("Failed to parse position: {:?}", parts[1]);
        }
        (0.0, 1, 0)
    }
}
