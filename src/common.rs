use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::time::Duration;
use clap::Parser;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    #[arg(short, long)]
    pub vcf: String,

    #[arg(short, long)]
    pub scoring: String,

    #[arg(long)]
    pub output: String,

    #[arg(long)]
    pub info: bool,
}

pub enum FileType {
    SingleSample,
    MultiSample,
    Gzipped,
}

impl FileType {
    pub fn detect(path: &str) -> io::Result<Self> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        let mut buffer = Vec::new();

        if path.ends_with(".gz") {
            Ok(FileType::Gzipped)
        } else {
            reader.read_until(b'\n', &mut buffer)?;
            if buffer.starts_with(b"##fileformat=VCF") {
                reader.read_until(b'\n', &mut buffer)?;
                while buffer.starts_with(b"##") {
                    buffer.clear();
                    reader.read_until(b'\n', &mut buffer)?;
                }
                let sample_count = buffer.split(|&b| b == b'\t').count() - 9;
                if sample_count > 1 {
                    Ok(FileType::MultiSample)
                } else {
                    Ok(FileType::SingleSample)
                }
            } else {
                Err(io::Error::new(io::ErrorKind::InvalidData, "Not a VCF file"))
            }
        }
    }
}

pub fn load_scoring_file(path: &str) -> io::Result<HashMap<(u8, u32), f32>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut effect_weights = HashMap::new();

    for line in reader.lines().filter(|l| l.as_ref().map_or(false, |s| !s.starts_with('#'))) {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 5 {
            if let (Ok(chr), Ok(pos), Ok(weight)) = (parts[0].parse::<u8>(), parts[1].parse::<u32>(), parts[4].parse::<f32>()) {
                effect_weights.insert((chr, pos), weight);
            }
        }
    }

    Ok(effect_weights)
}

pub fn output_results(args: &Args, score: f64, total_variants: usize, matched_variants: usize, duration: Duration, scoring_variants: usize) -> io::Result<()> {
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

pub fn print_info(score: f64, total_variants: usize, matched_variants: usize, scoring_variants: usize, duration: Duration) {
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
