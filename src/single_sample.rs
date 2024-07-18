use std::collections::HashMap;
use std::fs::File;
use rayon::prelude::*;
use memmap2::Mmap;
use std::str;
use std::io::{self, BufRead, BufReader};
use flate2::read::GzDecoder;
use noodles::vcf;
use noodles::bgzf;

pub fn calculate_polygenic_score(path: &str, effect_weights: &HashMap<(String, u32), f32>) -> io::Result<(f64, usize, usize)> {
    let mut reader = vcf::Reader::new(bgzf::Reader::new(File::open(path)?));

    let mut score = 0.0;
    let mut total_variants = 0;
    let mut matched_variants = 0;
    let mut debug_count = 0;

    // Read and discard the header
    reader.read_header()?;

    let mut record = vcf::Record::default();
    while reader.read_record(&mut record)? {
        total_variants += 1;

        if debug_count < 5 {
            println!("Raw VCF record: {:?}", record);
        }

        let chr = record.chromosome().to_string();
        let pos = record.position().get() as u32;

        debug_count += 1;
        if debug_count <= 5 {
            println!("Processing variant: chr={}, pos={}", chr, pos);
        }

        if let Some(&weight) = effect_weights.get(&(chr.clone(), pos)) {
            let genotype = record.genotypes().get("SAMPLE").and_then(|gt| gt.get(0));
            let allele_count = match genotype {
                Some(vcf::record::genotypes::sample::Value::Phased(0)) |
                Some(vcf::record::genotypes::sample::Value::Unphased(0)) => 0,
                Some(vcf::record::genotypes::sample::Value::Phased(1)) |
                Some(vcf::record::genotypes::sample::Value::Unphased(1)) => 1,
                Some(vcf::record::genotypes::sample::Value::Phased(2)) |
                Some(vcf::record::genotypes::sample::Value::Unphased(2)) => 2,
                _ => continue,
            };

            score += f64::from(weight) * allele_count as f64;
            matched_variants += 1;

            if debug_count <= 5 {
                println!("Matched variant: chr={}, pos={}, weight={}, allele_count={}", chr, pos, weight, allele_count);
            }
        }

        debug_count += 1;
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
