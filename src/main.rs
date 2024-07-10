use std::collections::HashMap;
use std::fs::File;
use std::time::Instant;
use rayon::prelude::*;
use clap::Parser;
use memmap2::Mmap;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(short, long)]
    vcf: String,

    #[arg(short, long)]
    scoring: String,

    #[arg(long)]
    output: String,

    #[arg(long)]
    info: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let start = Instant::now();

    let effect_weights = load_scoring_file(&args.scoring)?;
    let (score, total_variants, matched_variants) = calculate_polygenic_score(&args.vcf, &effect_weights)?;

    let duration = start.elapsed();

    output_results(&args, score, total_variants, matched_variants, duration, effect_weights.len())?;

    println!("Polygenic Score: {}", score);
    println!("Calculation time: {:?}", duration);

    if args.info {
        print_info(score, total_variants, matched_variants, effect_weights.len(), duration);
    }

    Ok(())
}

fn load_scoring_file(path: &str) -> Result<HashMap<(u8, u32), f32>, std::io::Error> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let mut effect_weights = HashMap::new();

    for line in mmap.split(|&b| b == b'\n').skip_while(|l| l.starts_with(b"#")) {
        let parts: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
        if parts.len() >= 5 {
            if let (Some(chr), Some(pos), Some(weight)) = (
                std::str::from_utf8(parts[0]).ok().and_then(|s| s.parse::<u8>().ok()),
                std::str::from_utf8(parts[1]).ok().and_then(|s| s.parse::<u32>().ok()),
                std::str::from_utf8(parts[4]).ok().and_then(|s| s.parse::<f32>().ok())
            ) {
                effect_weights.insert((chr, pos), weight);
            }
        }
    }

    Ok(effect_weights)
}

fn calculate_polygenic_score(path: &str, effect_weights: &HashMap<(u8, u32), f32>) -> Result<(f64, usize, usize), std::io::Error> {
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
        if let (Some(chr), Some(pos), Some(genotype)) = (parts.next(), parts.next(), parts.nth(7)) {
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

fn output_results(args: &Args, score: f64, total_variants: usize, matched_variants: usize, duration: std::time::Duration, scoring_variants: usize) -> std::io::Result<()> {
    let output = format!(
        "VCF_File\tScore_File\tPolygenic_Score\tCalculation_Time_Seconds\tTotal_Variants\tMatched_Variants\tScoring_Variants\n\
         {}\t{}\t{}\t{:.6}\t{}\t{}\t{}\n",
        args.vcf,
        args.scoring,
        score,
        duration.as_secs_f64(),
        total_variants,
        matched_variants,
        scoring_variants
    );

    std::fs::write(&args.output, output)
}

fn print_info(score: f64, total_variants: usize, matched_variants: usize, scoring_variants: usize, duration: std::time::Duration) {
    println!("\nDetailed Information:");
    println!("---------------------");
    println!("Total variants processed: {}", total_variants);
    println!("Variants in scoring file: {}", scoring_variants);
    println!("Matched variants: {}", matched_variants);
    println!("Match rate: {:.2}%", (matched_variants as f64 / scoring_variants as f64) * 100.0);
    println!("Polygenic Score: {}", score);
    println!("Calculation time: {:.6} seconds", duration.as_secs_f64());
    println!("Variants processed per second: {:.0}", total_variants as f64 / duration.as_secs_f64());
}