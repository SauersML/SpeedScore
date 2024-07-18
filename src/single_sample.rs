use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;

pub fn calculate_polygenic_score(path: &str, effect_weights: &HashMap<(String, u32), f32>) -> io::Result<(f64, usize, usize)> {
    let file = File::open(path)?;
    let reader = BufReader::with_capacity(1024 * 1024, MultiGzDecoder::new(file)); // 1MB buffer

    let lines: Vec<String> = reader.lines().collect::<io::Result<_>>()?;

    let (score, total_variants, matched_variants) = lines
        .par_iter()
        .enumerate()
        .filter(|(_, line)| !line.starts_with('#'))
        .map(|(index, line)| process_line(line, effect_weights, index))
        .reduce(
            || (0.0, 0, 0),
            |(score_a, total_a, matched_a), (score_b, total_b, matched_b)| {
                (score_a + score_b, total_a + total_b, matched_a + matched_b)
            },
        );

    Ok((score, total_variants, matched_variants))
}

fn process_line(line: &str, effect_weights: &HashMap<(String, u32), f32>, index: usize) -> (f64, usize, usize) {
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
    if let Ok(pos) = parts[1].parse::<u32>() {
        if index < 5 {
            println!("Processing variant: chr={}, pos={}", chr, pos);
        }

        if let Some(&weight) = effect_weights.get(&(chr.to_string(), pos)) {
            let genotype = parts[9];
            let allele_count = match genotype.chars().next() {
                Some('0') => 0,
                Some('1') => if genotype.chars().nth(2) == Some('1') { 2 } else { 1 },
                _ => return (0.0, 1, 0),
            };

            let score = f64::from(weight) * allele_count as f64;

            if index < 5 {
                println!("Matched variant: chr={}, pos={}, weight={}, allele_count={}", chr, pos, weight, allele_count);
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
