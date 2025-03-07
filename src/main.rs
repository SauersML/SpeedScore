use std::time::Instant;
use clap::Parser;
mod common;
mod single_sample;
mod multi_sample;
use common::{Args, FileType, load_scoring_file, output_results, print_info};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let start = Instant::now();
    let (effect_weights, scoring_chr_format) = load_scoring_file(&args.scoring)?;
    
    let file_type = FileType::detect(&args.vcf)?;
    
    let (score, total_variants, matched_variants, vcf_chr_format) = match file_type {
        FileType::SingleSample => {
            single_sample::calculate_polygenic_score(&args.vcf, &effect_weights)?
        },
        FileType::MultiSample => {
            let output_path = if args.output.is_empty() {
                format!("{}.csv", args.vcf)
            } else {
                args.output.clone()
            };
            multi_sample::calculate_polygenic_score_multi(
                &args.vcf,
                &effect_weights,
                &output_path,
                args.info
            )?
        },
    };

    let duration = start.elapsed();

    match file_type {
        FileType::SingleSample => {
            output_results(&args, score, total_variants, matched_variants, duration, effect_weights.len(), vcf_chr_format, scoring_chr_format)?;
            println!("Polygenic Score: {}", score);
        },
        FileType::MultiSample => {
            println!("Average Polygenic Score: {}", score);
        },
    }

    println!("Calculation time: {:?}", duration);
    println!("Total variants processed: {}", total_variants);
    println!("Matched variants: {}", matched_variants);

    if args.info {
        print_info(score, total_variants, matched_variants, effect_weights.len(), duration, vcf_chr_format, scoring_chr_format);
    }

    Ok(())
}
