use std::time::Instant;
use clap::Parser;
use std::path::Path;

mod common;
mod single_sample;
mod multi_sample;

use common::{Args, FileType, load_scoring_file, output_results, print_info};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let start = Instant::now();

    let effect_weights = load_scoring_file(&args.scoring)?;
    
    let file_type = FileType::detect(&args.vcf)?;
    
    let (score, total_variants, matched_variants) = match file_type {
        FileType::SingleSample => single_sample::calculate_polygenic_score(&args.vcf, &effect_weights)?,
        FileType::MultiSample | FileType::Gzipped => multi_sample::calculate_polygenic_score_multi(&args.vcf, &effect_weights)?,
    };

    let duration = start.elapsed();

    output_results(&args, score, total_variants, matched_variants, duration, effect_weights.len())?;

    println!("Polygenic Score: {}", score);
    println!("Calculation time: {:?}", duration);

    if args.info {
        print_info(score, total_variants, matched_variants, effect_weights.len(), duration);
    }

    Ok(())
}
