use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, BufRead, BufReader, Seek, SeekFrom};
use flate2::read::GzDecoder;

struct ResilienceReader<R: Read + Seek> {
    inner: R,
    recovery_attempts: usize,
}

impl<R: Read + Seek> ResilienceReader<R> {
    fn new(inner: R) -> Self {
        Self {
            inner,
            recovery_attempts: 0,
        }
    }
}

impl<R: Read + Seek> Read for ResilienceReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self.inner.read(buf) {
            Ok(0) if self.recovery_attempts < 3 => {
                self.recovery_attempts += 1;
                self.inner.seek(SeekFrom::Start(0))?;
                self.read(buf)
            }
            result => result,
        }
    }
}


pub fn calculate_polygenic_score_multi(path: &str, effect_weights: &HashMap<(u8, u32), f32>, debug: bool) -> io::Result<(f64, usize, usize)> {
    if debug {
        println!("Opening file: {}", path);
    }
    let file = File::open(path)?;
    
    let mut reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        if debug {
            println!("Detected gzipped file, using GzDecoder");
        }
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        if debug {
            println!("Using standard File reader");
        }
        Box::new(BufReader::new(file))
    };

    if debug {
        println!("Searching for header...");
    }
    let header = find_header(&mut reader)?;
    let sample_count = header.split('\t').count() - 9;

    if debug {
        println!("Header found: {}", header);
        println!("Sample count: {}", sample_count);
        println!("Attempting to read data after header:");
    }

    let mut line = String::new();
    let mut total_score = 0.0;
    let mut total_variants = 0;
    let mut total_matched = 0;

    while reader.read_line(&mut line)? > 0 {
        if !line.starts_with('#') {
            let (score, variants, matched) = process_line(&line, effect_weights, sample_count, debug);
            total_score += score;
            total_variants += variants;
            total_matched += matched;
            if debug && total_variants % 100_000 == 0 {
                println!("Processed {} variants", total_variants);
            }
        }
        line.clear();
    }

    if debug {
        println!("Total variants processed: {}", total_variants);
        if total_variants == 0 {
            println!("Warning: No variants were processed after the header!");
        }
    }

    Ok((total_score / sample_count as f64, total_variants, total_matched))
}

fn find_header(reader: &mut dyn BufRead) -> io::Result<String> {
    let mut line = String::new();
    loop {
        line.clear();
        reader.read_line(&mut line)?;
        if line.starts_with("#CHROM") {
            return Ok(line);
        }
    }
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
            if debug {
                println!("Variant matched. Chr: {}, Pos: {}, Score: {}, Matched: {}", chr, pos, score, matched);
            }
            return (score, 1, matched);
        } else if debug {
            println!("Variant not found in effect weights. Chr: {}, Pos: {}", chr, pos);
        }
    } else if debug {
        println!("Invalid chromosome or position: {}", line.trim());
    }
    (0.0, 1, 0)
}
