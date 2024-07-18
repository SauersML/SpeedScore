use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::time::Duration;
use clap::Parser;
use flate2::read::GzDecoder;

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
}

impl FileType {
    pub fn detect(path: &str) -> io::Result<Self> {
        let file = File::open(path)?;
        let mut reader: Box<dyn BufRead> = if path.ends_with(".gz") {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        let mut buffer = String::new();
        reader.read_line(&mut buffer)?;

        if !buffer.starts_with("##fileformat=VCF") {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "Not a VCF file"));
        }

        buffer.clear();
        while reader.read_line(&mut buffer)? > 0 {
            if buffer.starts_with("#CHROM") {
                let sample_count = buffer.split('\t').count() - 9;
                return Ok(if sample_count > 1 { FileType::MultiSample } else { FileType::SingleSample });
            }
            buffer.clear();
        }

        Err(io::Error::new(io::ErrorKind::InvalidData, "VCF header not found"))
    }
}



pub fn load_scoring_file(path: &str) -> io::Result<HashMap<(String, u32), f32>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut effect_weights = HashMap::new();
    let mut headers: Option<Vec<String>> = None;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        
        if headers.is_none() {
            headers = Some(line.split('\t').map(String::from).collect());
            continue;
        }

        let headers = headers.as_ref().unwrap();
        let parts: Vec<&str> = line.split('\t').collect();

        if parts.len() != headers.len() {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "Mismatch between header and data columns"));
        }

        let chr_index = headers.iter().position(|h| h == "chr_name").ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "Missing 'chr_name' column")
        })?;

        let pos_index = headers.iter().position(|h| h == "chr_position").ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "Missing 'chr_position' column")
        })?;

        let weight_index = headers.iter().position(|h| h == "effect_weight").ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "Missing 'effect_weight' column")
        })?;

        if let (Ok(chr), Ok(pos), Ok(weight)) = (
            parts[chr_index].parse::<String>(),
            parts[pos_index].parse::<u32>(),
            parts[weight_index].parse::<f32>()
        ) {
            effect_weights.insert((chr, pos), weight);
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
