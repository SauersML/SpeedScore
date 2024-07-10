use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader};
use flate2::read::GzDecoder;
use std::time::Instant;

const BUFFER_SIZE: usize = 1024 * 1024; // 1MB buffer

// Custom error type for our operations
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
            VcfError::InvalidFormat(msg) => write!(f, "Invalid VCF format: {}", msg),
            VcfError::Utf8Error(err) => write!(f, "UTF-8 error: {}", err),
        }
    }
}

impl std::error::Error for VcfError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            VcfError::Io(err) => Some(err),
            VcfError::InvalidFormat(_) => None,
            VcfError::Utf8Error(err) => Some(err),
        }
    }
}

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

// A trait for reading lines that works with both regular and gzipped files
trait LineReader {
    fn read_line(&mut self, buf: &mut String) -> Result<usize, VcfError>;
}

// Implement LineReader for BufReader<File>
impl LineReader for BufReader<File> {
    fn read_line(&mut self, buf: &mut String) -> Result<usize, VcfError> {
        buf.clear();
        Ok(io::BufRead::read_line(self, buf)?)
    }
}

// Implement LineReader for BufReader<GzDecoder<File>>
impl LineReader for BufReader<GzDecoder<File>> {
    fn read_line(&mut self, buf: &mut String) -> Result<usize, VcfError> {
        buf.clear();
        Ok(io::BufRead::read_line(self, buf)?)
    }
}

// A struct to hold our VCF reader and its metadata
struct VcfReader {
    reader: Box<dyn LineReader>,
    sample_count: usize,
}

impl VcfReader {
    fn new(path: &str) -> Result<Self, VcfError> {
        let file = File::open(path)?;
        let reader: Box<dyn LineReader> = if path.ends_with(".gz") {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        let mut vcf_reader = VcfReader { reader, sample_count: 0 };
        vcf_reader.find_header()?;
        Ok(vcf_reader)
    }

    fn find_header(&mut self) -> Result<(), VcfError> {
        let mut line = String::new();
        loop {
            self.reader.read_line(&mut line)?;
            if line.starts_with("#CHROM") {
                self.sample_count = line.split_whitespace().count() - 9;
                return Ok(());
            } else if !line.starts_with('#') {
                return Err(VcfError::InvalidFormat("VCF header not found".to_string()));
            }
        }
    }
}

pub fn calculate_polygenic_score_multi(
    path: &str,
    effect_weights: &HashMap<(u8, u32), f32>,
    debug: bool
) -> Result<(f64, usize, usize), VcfError> {
    let start_time = Instant::now();
    
    if debug {
        println!("Opening file: {}", path);
    }

    let mut vcf_reader = VcfReader::new(path)?;

    if debug {
        println!("VCF data start found.");
        println!("Sample count: {}", vcf_reader.sample_count);
        println!("Processing variants...");
    }

    let mut line = String::new();
    let mut total_score = 0.0;
    let mut total_variants = 0;
    let mut total_matched = 0;

    while vcf_reader.reader.read_line(&mut line)? > 0 {
        if !line.starts_with('#') {
            let (score, variants, matched) = process_line(&line, effect_weights, vcf_reader.sample_count, debug);
            total_score += score;
            total_variants += variants;
            total_matched += matched;

            if debug && total_variants % 100_000 == 0 {
                println!("Processed {} variants", total_variants);
            }
        }
    }

    let duration = start_time.elapsed();

    if debug {
        println!("Total variants processed: {}", total_variants);
        println!("Matched variants: {}", total_matched);
        println!("Processing time: {:?}", duration);
        if total_variants == 0 {
            println!("Warning: No variants were processed!");
        }
    }

    Ok((total_score / vcf_reader.sample_count as f64, total_variants, total_matched))
}

fn process_line(line: &str, effect_weights: &HashMap<(u8, u32), f32>, sample_count: usize, debug: bool) -> (f64, usize, usize) {
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
            
            if debug && matched > 0 {
                println!("Matched variant at Chr {}, Pos {}. Score: {}, Matched: {}", chr, pos, score, matched);
            }
            
            return (score, 1, matched);
        }
    }
    
    (0.0, 1, 0)
}
