use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::time::Duration;
use clap::Parser;
use flate2::read::GzDecoder;
use std::cell::RefCell;
use std::rc::Rc;

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

pub struct ChromosomeFormat {
    pub has_chr_prefix: bool,
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


pub fn load_scoring_file(
    path: &str
) -> io::Result<(HashMap<(String, u32), (String, f32)>, bool)> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut effect_weights: HashMap<(String, u32), (String, f32)> = HashMap::new();
    let mut headers: Option<Vec<String>> = None;
    let mut scoring_chr_format = false;

    let mut count = 0;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }

        // First non‐comment line is assumed to be headers
        if headers.is_none() {
            headers = Some(line.split('\t').map(String::from).collect());
            continue;
        }

        let headers = headers.as_ref().unwrap();
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() != headers.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Mismatch between header and data columns"
            ));
        }

        // Find column indices for chr, position, effect_allele, effect_weight
        let chr_index = headers.iter().position(|h| h == "chr_name").ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "Missing 'chr_name' column")
        })?;

        let pos_index = headers.iter().position(|h| h == "chr_position").ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "Missing 'chr_position' column")
        })?;

        let allele_index = headers.iter().position(|h| h == "effect_allele").ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "Missing 'effect_allele' column")
        })?;

        let weight_index = headers.iter().position(|h| h == "effect_weight").ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "Missing 'effect_weight' column")
        })?;

        let chr = parts[chr_index].to_string();
        let pos = parts[pos_index].parse::<u32>().map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, "Invalid numeric position")
        })?;
        let allele = parts[allele_index].to_string();  // e.g., "A", "T", etc.
        let weight = parts[weight_index].parse::<f32>().map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, "Invalid numeric weight")
        })?;

        // Check if our first line uses 'chr' prefix
        if count == 0 {
            scoring_chr_format = chr.starts_with("chr");
        }

        // Normalize chromosome (remove leading "chr")
        let normalized_chr = chr.trim_start_matches("chr").to_string();

        // Store (effect_allele, effect_weight)
        effect_weights.insert((normalized_chr, pos), (allele, weight));
        count += 1;

        if count <= 5 {
            println!(
                "Loaded scoring data example: chr={}, pos={}, allele={}, weight={}",
                chr, pos, allele, weight
            );
        }
    }

    println!("Total scoring entries loaded: {}", effect_weights.len());
    Ok((effect_weights, scoring_chr_format))
}


pub fn output_results(args: &Args, score: f64, total_variants: usize, matched_variants: usize, duration: Duration, scoring_variants: usize, vcf_chr_format: bool, scoring_chr_format: bool) -> io::Result<()> {
    let output = format!(
        "VCF_File\tScore_File\tPolygenic_Score\tCalculation_Time_Seconds\tTotal_Variants\tMatched_Variants\tScoring_Variants\tVCF_Chr_Format\tScoring_Chr_Format\n\
         {}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{}\n",
        args.vcf,
        args.scoring,
        score,
        duration.as_secs_f64(),
        total_variants,
        matched_variants,
        scoring_variants,
        vcf_chr_format,
        scoring_chr_format
    );

    std::fs::write(&args.output, output)
}

pub fn print_info(score: f64, total_variants: usize, matched_variants: usize, scoring_variants: usize, duration: Duration, vcf_chr_format: bool, scoring_chr_format: bool) {
    println!("\nDetailed Information:");
    println!("---------------------");
    println!("Total variants processed: {}", total_variants);
    println!("Variants in scoring file: {}", scoring_variants);
    println!("Matched variants: {}", matched_variants);
    println!("Match rate: {:.2}%", (matched_variants as f64 / scoring_variants as f64) * 100.0);
    println!("Polygenic Score: {}", score);
    println!("Calculation time: {:.6} seconds", duration.as_secs_f64());
    println!("Variants processed per second: {:.0}", total_variants as f64 / duration.as_secs_f64());
    println!("VCF chromosome format: {}", if vcf_chr_format { "chr" } else { "no chr" });
    println!("Scoring file chromosome format: {}", if scoring_chr_format { "chr" } else { "no chr" });
}
