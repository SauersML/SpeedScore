use std::collections::HashMap;
use std::fs::File;
use rayon::prelude::*;
use memmap2::Mmap;
use std::io::{self, BufRead};
use noodles::vcf;

pub fn calculate_polygenic_score(path: &str, effect_weights: &HashMap<(u8, u32), f32>) -> io::Result<(f64, usize, usize)> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let chunk_size = 1024 * 1024 * 10; // 10MB chunks
    let num_chunks = (mmap.len() + chunk_size - 1) / chunk_size;
    let results: Vec<_> = (0..num_chunks)
        .into_par_iter()
        .map(|i| {
            let start = i * chunk_size;
            let end = (start + chunk_size).min(mmap.len());
            process_chunk(&mmap[start..end], effect_weights)
        })
        .collect();
    let (score, total_variants, matched_variants) = results.into_iter()
        .fold((0.0, 0, 0), |acc, x| (acc.0 + x.0, acc.1 + x.1, acc.2 + x.2));
    Ok((score, total_variants, matched_variants))
}

fn process_chunk(chunk: &[u8], effect_weights: &HashMap<(u8, u32), f32>) -> (f64, usize, usize) {
    let mut score = 0.0;
    let mut total_variants = 0;
    let mut matched_variants = 0;
    for line in chunk.split(|&b| b == b'\n') {
        if line.is_empty() || line[0] == b'#' {
            continue;
        }
        total_variants += 1;
        let mut parts = line.split(|&b| b == b'\t');
        if let (Some(chr), Some(pos), Some(_ref), Some(_alt), Some(genotype)) = (parts.next(), parts.next(), parts.next(), parts.next(), parts.nth(5)) {
            if let (Some(chr), Some(pos)) = (
                std::str::from_utf8(chr).ok().and_then(|s| s.parse::<u8>().ok()),
                std::str::from_utf8(pos).ok().and_then(|s| s.parse::<u32>().ok())
            ) {
                if let Some(&weight) = effect_weights.get(&(chr, pos)) {
                    let allele_count = match genotype.get(0) {
                        Some(b'0') => 0,
                        Some(b'1') => if genotype.get(2) == Some(&b'1') { 2 } else { 1 },
                        _ => continue,
                    };
                    score += f64::from(weight) * allele_count as f64;
                    matched_variants += 1;
                }
            }
        }
    }
    (score, total_variants, matched_variants)
}

pub fn debug_first_lines(path: &str, num_lines: usize) -> io::Result<()> {
    let mut reader = vcf::Reader::new(File::open(path)?);

    println!("VCF header:");
    match reader.read_header() {
        Ok(header) => println!("{}", header),
        Err(e) => println!("Error reading header: {:?}", e),
    }

    println!("\nFirst {} data lines:", num_lines);
    for (i, result) in reader.records().take(num_lines).enumerate() {
        match result {
            Ok(record) => println!("Line {}: {:?}", i, record),
            Err(e) => println!("Error reading line {}: {:?}", i, e),
        }
    }

    Ok(())
}
