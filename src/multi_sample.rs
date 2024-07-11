use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, BufRead, BufReader};
use std::time::Instant;
use flate2::read::MultiGzDecoder;
use indicatif::{ProgressBar, ProgressStyle, MultiProgress};




#[derive(Debug)]
pub enum VcfError {
    Io(io::Error),
    InvalidFormat(String),
    Utf8Error(std::string::FromUtf8Error),
}

impl std::fmt::Display for VcfError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            VcfError::Io(err) => write!(f, "I/O error: {}", err),
            VcfError::InvalidFormat(msg) => write!(f, "Invalid format: {}", msg),
            VcfError::Utf8Error(err) => write!(f, "UTF-8 error: {}", err),
        }
    }
}

impl std::error::Error for VcfError {}

impl From<io::Error> for VcfError {
    fn from(error: io::Error) -> Self {
        VcfError::Io(error)
    }
}

impl From<std::string::FromUtf8Error> for VcfError {
    fn from(error: std::string::FromUtf8Error) -> Self {
        VcfError::Utf8Error(error)
    }
}

struct VcfReader<R: Read> {
    reader: BufReader<R>,
    sample_count: usize,
}

impl<R: Read> VcfReader<R> {
    fn new(reader: R) -> Result<Self, VcfError> {
        let mut vcf_reader = VcfReader {
            reader: BufReader::new(reader),
            sample_count: 0,
        };
        vcf_reader.find_header()?;
        Ok(vcf_reader)
    }

    fn find_header(&mut self) -> Result<(), VcfError> {
        let mut line = String::new();
        loop {
            if self.reader.read_line(&mut line)? == 0 {
                return Err(VcfError::InvalidFormat("VCF header not found".to_string()));
            }
            if line.starts_with("#CHROM") {
                self.sample_count = line.split_whitespace().count() - 9;
                return Ok(());
            }
            line.clear();
        }
    }

    fn read_line(&mut self, buf: &mut String) -> Result<usize, VcfError> {
        buf.clear();
        Ok(self.reader.read_line(buf)?)
    }
}

fn open_vcf_reader(path: &str) -> Result<VcfReader<MultiGzDecoder<File>>, VcfError> {
    let file = File::open(path)?;
    let decoder = MultiGzDecoder::new(file);
    VcfReader::new(decoder)
}













pub fn calculate_polygenic_score_multi(
    path: &str,
    effect_weights: &HashMap<(u8, u32), f32>,
    debug: bool
) -> Result<(f64, usize, usize), VcfError> {
    let start_time = Instant::now();

    println!("Opening file: {}", path);
    println!("Effect weights loaded: {} variants", effect_weights.len());

    let mut vcf_reader = open_vcf_reader(path)?;

    println!("VCF data start found.");
    println!("Sample count: {}", vcf_reader.sample_count);
    println!("Processing variants...");

    // Set up progress bars
    let multi_progress = MultiProgress::new();
    let row_pb = multi_progress.add(ProgressBar::new_spinner());
    row_pb.set_style(ProgressStyle::default_spinner()
        .template("{spinner:.green} [{elapsed_precise}] {msg}")
        .unwrap());
    row_pb.set_message("Starting processing...");

    let sample_pb = multi_progress.add(ProgressBar::new(vcf_reader.sample_count as u64));
    sample_pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} samples processed {msg}")
        .unwrap());

    let mut line = String::new();
    let mut total_score = 0.0;
    let mut total_variants = 0;
    let mut total_matched = 0;
    let mut in_data_section = false;
    let mut lines_processed = 0;

    while vcf_reader.read_line(&mut line)? > 0 {
        lines_processed += 1;

        if line.starts_with("#CHROM") {
            in_data_section = true;
            row_pb.set_message("Found #CHROM line, processing data...");
            continue;
        }

        if in_data_section {
            if lines_processed <= 1000 && lines_processed % 100 == 0 {
                let truncated_line: String = line.split('\t').take(14).collect::<Vec<_>>().join("\t");
                println!("\nLine {}: {}", lines_processed, truncated_line);
            }

            let (score, variants, matched) = process_line(&line, effect_weights, vcf_reader.sample_count, debug);
            total_score += score;
            total_variants += variants;
            total_matched += matched;

            if lines_processed % 1000 == 0 {
                row_pb.set_message(format!("Processed {} rows, {} variants, matched {}", lines_processed, total_variants, total_matched));
            }

            if lines_processed % 10_000 == 0 {
                let samples_processed = (lines_processed / 10_000).min(vcf_reader.sample_count);
                sample_pb.set_position(samples_processed as u64);
                sample_pb.set_message(format!("Avg score: {:.6}", total_score / vcf_reader.sample_count as f64));
            }

            if lines_processed % 100_000 == 0 {
                let avg_score_per_sample = total_score / vcf_reader.sample_count as f64;
                println!("\nProgress update:");
                println!("Lines processed: {}", lines_processed);
                println!("Variants processed: {}", total_variants);
                println!("Matched variants: {}", total_matched);
                println!("Current average score per sample: {:.6}", avg_score_per_sample);
            }
        }
    }

    row_pb.finish_with_message("Processing complete");
    sample_pb.finish_with_message("All samples processed");

    let duration = start_time.elapsed();

    println!("\nFinished reading all lines.");
    println!("Total lines processed: {}", lines_processed);
    println!("Total variants processed: {}", total_variants);
    println!("Matched variants: {}", total_matched);
    println!("Final average score per sample: {:.6}", total_score / vcf_reader.sample_count as f64);
    println!("Processing time: {:?}", duration);

    if total_variants == 0 {
        println!("Warning: No variants were processed!");
    }

    Ok((total_score / vcf_reader.sample_count as f64, total_variants, total_matched))
}










fn process_line(line: &str, effect_weights: &HashMap<(u8, u32), f32>, sample_count: usize, debug: bool) -> (f64, usize, usize) {
    let mut parts = line.split('\t');
    let chr = parts.next().and_then(|s| s.parse::<u8>().ok());
    let pos = parts.next().and_then(|s| s.parse::<u32>().ok());

    if debug {
        println!("Parsed chr: {:?}, pos: {:?}", chr, pos);
    }

    if let (Some(chr), Some(pos)) = (chr, pos) {
        if let Some(&weight) = effect_weights.get(&(chr, pos)) {
            if debug {
                println!("Found matching weight for chr: {}, pos: {}", chr, pos);
            }
            let genotypes: Vec<&str> = parts.skip(7).take(sample_count).collect();
            if debug {
                println!("Genotypes: {:?}", genotypes.iter().take(5).collect::<Vec<_>>());
            }
            let (score, matched) = genotypes.iter()
                .map(|&gt| match gt.chars().next() {
                    Some('0') => (0.0, 1),
                    Some('1') => (f64::from(weight), 1),
                    _ => (0.0, 0),
                })
                .fold((0.0, 0), |acc, x| (acc.0 + x.0, acc.1 + x.1));
            
            if debug && matched > 0 {
                println!("Matched variant at Chr {}, Pos {}. Score: {}, Matched: {}", chr, pos, score, matched);
            }
            
            return (score, 1, matched);
        } else {
            if debug {
                println!("No matching weight for chr: {}, pos: {}", chr, pos);
            }
        }
    } else {
        if debug {
            println!("Could not parse chr or pos from line: {}", line);
        }
    }
    
    (0.0, 1, 0)
}
