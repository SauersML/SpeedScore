use std::time::Instant;
use clap::Parser;
mod common;
mod single_sample;
mod multi_sample;
use common::{Args, FileType, load_scoring_file, output_results, print_info};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let start = Instant::now();
    let effect_weights = load_scoring_file(&args.scoring)?;
    
    println!("Loaded {} variants from scoring file", effect_weights.len());
    println!("First 5 entries in effect_weights:");
    for (i, ((chr, pos), weight)) in effect_weights.iter().take(5).enumerate() {
        println!("  {}: Chr {}, Pos {}, Weight {}", i+1, chr, pos, weight);
    }

    println!("\nFirst 5 lines of VCF file:");
    single_sample::debug_first_lines(&args.vcf, 5)?;

    let file_type = FileType::detect(&args.vcf)?;
    
    let (score, total_variants, matched_variants) = match file_type {
        FileType::SingleSample => {
            let (score, total, matched) = single_sample::calculate_polygenic_score(&args.vcf, &effect_weights)?;
            (score, total, matched)
        },
        FileType::MultiSample => {
            let (score, total, matched) = single_sample::calculate_polygenic_score(&args.vcf, &effect_weights)?;
            (score, total, matched)
        },
        FileType::MultiSample => {
            let output_path = if args.output.is_empty() {
                format!("{}.csv", args.vcf)
            } else {
                args.output.clone()
            };
            let (avg_score, total, matched) = multi_sample::calculate_polygenic_score_multi(
                &args.vcf,
                &effect_weights,
                &output_path,
                args.info
            )?;
            println!("Multi-sample results written to: {}", output_path);
            (avg_score, total, matched)
        },
    };

    let duration = start.elapsed();

    match file_type {
        FileType::SingleSample => {
            output_results(&args, score, total_variants, matched_variants, duration, effect_weights.len())?;
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
        print_info(score, total_variants, matched_variants, effect_weights.len(), duration);
    }

    Ok(())
}
