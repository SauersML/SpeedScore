use std::collections::HashMap;
use std::fs::File;
use std::io;
use rayon::prelude::*;
use memmap2::Mmap;

pub fn calculate_polygenic_score(path: &str, effect_weights: &HashMap<(String, u32), f32>) -> io::Result<(f64, usize, usize)> {
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
            if let (Ok(chr), Ok(pos)) = (
                std::str::from_utf8(chr).map(|s| s.trim().to_string()),
                std::str::from_utf8(pos).and_then(|s| Ok(s.trim().parse::<u32>().ok()))
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
