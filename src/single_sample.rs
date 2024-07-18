use std::collections::HashMap;
use std::fs::File;
use rayon::prelude::*;
use memmap2::Mmap;
use std::str;
use noodles::bgzf;
use noodles::vcf::{self, record::GenotypeField};
use std::collections::HashMap;
use std::io::{self, BufRead, BufReader};
use flate2::read::MultiGzDecoder;

pub fn calculate_polygenic_score(path: &str, effect_weights: &HashMap<(String, u32), f32>) -> io::Result<(f64, usize, usize)> {
    let file = File::open(path)?;
    let reader = BufReader::new(MultiGzDecoder::new(file));

    let mut score = 0.0;
    let mut total_variants = 0;
    let mut matched_variants = 0;

    for (index, line) in reader.lines().enumerate() {
        let line = line?;
        
        if line.starts_with('#') {
            continue;
        }

        total_variants += 1;

        if index < 5 {
            println!("Raw VCF line {}: {}", index + 1, line);
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 10 {
            if index < 5 {
                println!("Invalid VCF line format: {:?}", parts);
            }
            continue;
        }

        let chr = parts[0].to_string();
        if let Ok(pos) = parts[1].parse::<u32>() {
            if index < 5 {
                println!("Processing variant: chr={}, pos={}", chr, pos);
            }

            if let Some(&weight) = effect_weights.get(&(chr.clone(), pos)) {
                let genotype = parts[9];
                let allele_count = match genotype.chars().next() {
                    Some('0') => 0,
                    Some('1') => if genotype.chars().nth(2) == Some('1') { 2 } else { 1 },
                    _ => continue,
                };

                score += f64::from(weight) * allele_count as f64;
                matched_variants += 1;

                if index < 5 {
                    println!("Matched variant: chr={}, pos={}, weight={}, allele_count={}", chr, pos, weight, allele_count);
                }
            }
        } else {
            if index < 5 {
                println!("Failed to parse position: {:?}", parts[1]);
            }
        }
    }

    Ok((score, total_variants, matched_variants))
}







fn process_chunk(chunk: &[u8], effect_weights: &HashMap<(String, u32), f32>) -> (f64, usize, usize) {
    let mut score = 0.0;
    let mut total_variants = 0;
    let mut matched_variants = 0;

    let mut debug_count = 0;
    for line in chunk.split(|&b| b == b'\n') {
        if line.is_empty() || line[0] == b'#' {
            continue;
        }
        total_variants += 1;

        if debug_count < 5 {
            println!("Raw VCF line: {}", String::from_utf8_lossy(line));
        }

        let mut parts = line.split(|&b| b == b'\t');
        if let (Some(chr), Some(pos), Some(genotype)) = (parts.next(), parts.next(), parts.nth(7)) {
            if let (Some(chr), Some(pos)) = (
                str::from_utf8(chr).ok().map(|s| s.trim().to_string()),
                str::from_utf8(pos).ok().and_then(|s| s.trim().parse::<u32>().ok())
            ) {
                debug_count += 1;
                if debug_count <= 5 {
                    println!("Processing variant: chr={}, pos={}", chr, pos);
                }
                if let Some(&weight) = effect_weights.get(&(chr.clone(), pos)) {
                    let allele_count = match genotype.get(0) {
                        Some(b'0') => 0,
                        Some(b'1') => if genotype.get(2) == Some(&b'1') { 2 } else { 1 },
                        _ => continue,
                    };
                    score += f64::from(weight) * allele_count as f64;
                    matched_variants += 1;
                    if debug_count <= 5 {
                        println!("Matched variant: chr={}, pos={}, weight={}, allele_count={}", chr, pos, weight, allele_count);
                    }
                }
            } else {
                if debug_count < 5 {
                    println!("Failed to parse chr or pos: chr={:?}, pos={:?}", 
                             String::from_utf8_lossy(chr), String::from_utf8_lossy(pos));
                }
            }
        } else {
            if debug_count < 5 {
                println!("Failed to extract chr, pos, or genotype from VCF line");
            }
        }

        debug_count += 1;
    }

    (score, total_variants, matched_variants)
}
