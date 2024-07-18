use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{self, Read, BufRead, BufReader, Write};
use std::time::Instant;
use std::path::Path;
use flate2::read::MultiGzDecoder;
use indicatif::{ProgressBar, ProgressStyle};

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
    sample_names: Vec<String>,
}

impl<R: Read> VcfReader<R> {
    fn new(reader: R) -> Result<Self, VcfError> {
        let mut vcf_reader = VcfReader {
            reader: BufReader::new(reader),
            sample_names: Vec::new(),
        };
        vcf_reader.find_header()?;
        Ok(vcf_reader)
    }

    fn find_header(&mut self) -> Result<(), VcfError> {
        for line in self.reader.by_ref().lines() {
            let line = line?;
            if line.starts_with("#CHROM") {
                self.sample_names = line.split_whitespace().skip(9).map(String::from).collect();
                return Ok(());
            }
        }
        Err(VcfError::InvalidFormat("VCF header not found".to_string()))
    }

    fn read_line(&mut self, buf: &mut String) -> Result<usize, VcfError> {
        buf.clear();
        Ok(self.reader.read_line(buf)?)
    }
}

#[derive(Clone, Default)]
struct SampleData {
    score: f64,
    matched_variants: usize,
    total_variants: usize,
}




fn open_vcf_reader(path: &str) -> Result<BufReader<MultiGzDecoder<File>>, VcfError> {
    let file = File::open(path).map_err(VcfError::Io)?;
    let decoder = MultiGzDecoder::new(file);
    Ok(BufReader::with_capacity(1024 * 1024, decoder)) // 1MB buffer
}

pub fn calculate_polygenic_score_multi(
    vcf_path: &str,
    effect_weights: &HashMap<(u8, u32), f32>,
    output_path: &str,
    debug: bool
) -> Result<(f64, usize, usize), VcfError> {
    let start_time = Instant::now();

    println!("Opening file: {}", vcf_path);
    println!("Effect weights loaded: {} variants", effect_weights.len());

    let mut reader = open_vcf_reader(vcf_path)?;
    let mut header_line = String::new();
    let mut sample_names = Vec::new();

    // Find the header
    loop {
        reader.read_line(&mut header_line)?;
        if header_line.starts_with("#CHROM") {
            sample_names = header_line.split_whitespace().skip(9).map(String::from).collect();
            break;
        }
        header_line.clear();
    }

    println!("VCF data start found.");
    println!("Sample count: {}", sample_names.len());
    println!("Processing variants...");

    let pb = ProgressBar::new_spinner();
    pb.set_style(ProgressStyle::default_spinner()
        .template("{spinner:.green} [{elapsed_precise}] {msg}")
        .unwrap());
    pb.set_message("Processing...");

    let mut buffer = Vec::new();
    let mut sample_data: Vec<SampleData> = vec![SampleData::default(); sample_names.len()];
    let mut lines_processed = 0;
    let mut last_chr = 0;
    let mut last_pos = 0;

    loop {
        buffer.clear();
        let num_lines = reader.read_until(b'\n', &mut buffer)?;
        if num_lines == 0 {
            break;
        }

        lines_processed += 1;

        if !buffer.starts_with(&[b'#']) {
            let result = process_chunk(&buffer, effect_weights, &mut sample_data, debug);
            if let Some((chr, pos)) = result {
                if debug && (chr != last_chr || pos > last_pos + 20_000_000) {
                    pb.suspend(|| {
                        println!("\rProcessed up to Chr {}, Pos {:.2}M", chr, pos as f64 / 1_000_000.0);
                        io::stdout().flush().unwrap();
                    });
                    last_chr = chr;
                    last_pos = pos;
                }
            }
        }

        if lines_processed % 100_000 == 0 {
            let lines_in_k = lines_processed / 1000;
            let variants = sample_data.iter().map(|sd| sd.total_variants).sum::<usize>();
            let matched = sample_data.iter().map(|sd| sd.matched_variants).sum::<usize>();
            pb.set_message(format!(
                "{}K lines, {}K variants, {}K matched",
                lines_in_k,
                variants / 1000,
                matched / 1000
            ));
        }
    }

    pb.finish_with_message("Processing complete");

    let duration = start_time.elapsed();

    write_csv_output(output_path, vcf_path, &sample_names, &sample_data, duration)?;

    let avg_score = sample_data.iter().map(|sd| sd.score).sum::<f64>() / sample_data.len() as f64;
    let total_variants = sample_data.iter().map(|sd| sd.total_variants).sum();
    let matched_variants = sample_data.iter().map(|sd| sd.matched_variants).sum();

    println!("\nFinished processing.");
    println!("Total lines processed: {:.3}K", lines_processed as f64 / 1000.0);
    println!("Results written to: {}", output_path);
    println!("Processing time: {:?}", duration);

    Ok((avg_score, total_variants, matched_variants))
}

fn process_chunk(chunk: &[u8], effect_weights: &HashMap<(u8, u32), f32>, sample_data: &mut [SampleData], _debug: bool) -> Option<(u8, u32)> {
    let mut last_chr = 0;
    let mut last_pos = 0;

    for line in chunk.split(|&b| b == b'\n') {
        if line.is_empty() {
            continue;
        }

        let mut parts = line.split(|&b| b == b'\t');
        let chr = parts.next().and_then(|s| std::str::from_utf8(s).ok()?.parse::<u8>().ok()).unwrap_or(0);
        let pos = parts.next().and_then(|s| std::str::from_utf8(s).ok()?.parse::<u32>().ok()).unwrap_or(0);

        if let Some(&weight) = effect_weights.get(&(chr, pos)) {
            let genotypes = parts.skip(7);
            for (sample, genotype) in sample_data.iter_mut().zip(genotypes) {
                sample.total_variants += 1;
                match genotype.first() {
                    Some(b'0') => sample.matched_variants += 1,
                    Some(b'1') => {
                        sample.score += f64::from(weight);
                        sample.matched_variants += 1;
                    }
                    _ => {}
                }
            }
        }

        last_chr = chr;
        last_pos = pos;
    }

    Some((last_chr, last_pos))
}

fn write_csv_output(
    output_path: &str,
    vcf_path: &str,
    sample_names: &[String],
    sample_data: &[SampleData],
    duration: std::time::Duration
) -> Result<(), VcfError> {
    let path = Path::new(output_path);
    let prefix = path.parent().unwrap_or_else(|| Path::new(""));
    std::fs::create_dir_all(prefix).map_err(VcfError::Io)?;

    let mut file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(output_path)
        .map_err(VcfError::Io)?;

    writeln!(file, "VCF_File,Sample_Name,Polygenic_Score,Calculation_Time_Seconds,Total_Variants,Matched_Variants")
        .map_err(VcfError::Io)?;

    for (name, data) in sample_names.iter().zip(sample_data.iter()) {
        writeln!(
            file,
            "{},{},{:.6},{:.6},{},{}",
            vcf_path,
            name,
            data.score,
            duration.as_secs_f64(),
            data.total_variants,
            data.matched_variants
        ).map_err(VcfError::Io)?;
    }

    Ok(())
}
