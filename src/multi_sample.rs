use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use rayon::prelude::*;
use flate2::read::GzDecoder;

pub fn calculate_polygenic_score_multi(path: &str, effect_weights: &HashMap<(u8, u32), f32>) -> io::Result<(f64, usize, usize)> {
    let file = File::open(path)?;
    let reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut lines = reader.lines();
    let header = find_header(&mut lines)?;
    let sample_count = header.split('\t').count() - 9;

    let (total_score, total_variants, total_matched) = lines
        .collect::<io::Result<Vec<_>>>()?
        .into_par_iter()
        .map(|line| process_line(&line, effect_weights, sample_count))
        .reduce(|| (0.0, 0, 0),
                |acc, x| (acc.0 + x.0, acc.1 + x.1, acc.2 + x.2));

    Ok((total_score / sample_count as f64, total_variants, total_matched))
}

fn find_header<B: BufRead>(lines: &mut std::io::Lines<B>) -> io::Result<String> {
    for line in lines {
        let line = line?;
        if line.starts_with("#CHROM") {
            return Ok(line);
        }
    }
    Err(io::Error::new(io::ErrorKind::InvalidData, "Header not found"))
}

fn process_line(line: &str, effect_weights: &HashMap<(u8, u32), f32>, sample_count: usize) -> (f64, usize, usize) {
    let mut parts = line.split('\t');
    let chr = parts.next().and_then(|s| s.parse::<u8>().ok());
    let pos = parts.next().and_then(|s| s.parse::<u32>().ok());

    if let (Some(chr), Some(pos)) = (chr, pos) {
        if let Some(&weight) = effect_weights.get(&(chr, pos)) {
            let genotypes: Vec<&str> = parts.skip(7).take(sample_count).collect();
            let (score, matched) = genotypes.iter()
                .map(|&gt| match gt.chars().next() {
                    Some('0') => (0.0, 1),
                    Some('1') => (f64::from(weight), 1),
                    _ => (0.0, 0),
                })
                .fold((0.0, 0), |acc, x| (acc.0 + x.0, acc.1 + x.1));
            return (score, 1, matched);
        }
    }
    (0.0, 1, 0)
}
