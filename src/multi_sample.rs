use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, BufReader};
use flate2::read::GzDecoder;

pub fn calculate_polygenic_score_multi(path: &str, effect_weights: &HashMap<(u8, u32), f32>, debug: bool) -> io::Result<(f64, usize, usize)> {
    if debug {
        println!("Opening file: {}", path);
    }
    let file = File::open(path)?;
    
    let mut reader: Box<dyn Read> = if path.ends_with(".gz") {
        if debug {
            println!("Detected gzipped file, using GzDecoder");
        }
        Box::new(GzDecoder::new(file))
    } else {
        if debug {
            println!("Using standard File reader");
        }
        Box::new(file)
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

    let mut buffer = [0; 1024 * 1024]; // 1MB buffer
    let mut line_buffer = Vec::new();
    let mut line_count = 0;
    let mut total_score = 0.0;
    let mut total_variants = 0;
    let mut total_matched = 0;

    loop {
        match reader.read(&mut buffer) {
            Ok(0) => break, // End of file
            Ok(n) => {
                if debug && line_count == 0 {
                    println!("Read {} bytes", n);
                }
                for &byte in &buffer[..n] {
                    if byte == b'\n' {
                        if !line_buffer.is_empty() {
                            let line = String::from_utf8_lossy(&line_buffer);
                            if debug && line_count < 5 {
                                println!("Line {}: {}", line_count + 1, line);
                            }
                            let (score, variants, matched) = process_line(&line, effect_weights, sample_count, debug);
                            total_score += score;
                            total_variants += variants;
                            total_matched += matched;
                            line_count += 1;
                            if debug && line_count % 100_000 == 0 {
                                println!("Processed {} lines", line_count);
                            }
                        }
                        line_buffer.clear();
                    } else {
                        line_buffer.push(byte);
                    }
                }
            }
            Err(e) => {
                if debug {
                    println!("Error reading file: {}", e);
                }
                return Err(e);
            }
        }
    }

    if debug {
        println!("Total lines processed: {}", line_count);
        if line_count == 0 {
            println!("Warning: No lines were processed after the header!");
        }
    }

    Ok((total_score / sample_count as f64, total_variants, total_matched))
}

fn find_header(reader: &mut dyn Read) -> io::Result<String> {
    let mut buffer = Vec::new();
    let mut byte = [0u8; 1];
    
    loop {
        reader.read_exact(&mut byte)?;
        buffer.push(byte[0]);
        
        if byte[0] == b'\n' {
            let line = String::from_utf8_lossy(&buffer);
            if line.starts_with("#CHROM") {
                return Ok(line.into_owned());
            }
            buffer.clear();
        }
    }
}

fn process_line(line: &str, effect_weights: &HashMap<(u8, u32), f32>, sample_count: usize, debug: bool) -> (f64, usize, usize) {
    let mut parts = line.split('\t');
    let chr = parts.next().and_then(|s| s.parse::<u8>().ok());
    let pos = parts.next().and_then(|s| s.parse::<u32>().ok());

    if let (Some(chr), Some(pos)) = (chr, pos) {
        if debug {
            println!("Processing variant: Chr {}, Pos {}", chr, pos);
        }
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
                println!("Variant matched. Score: {}, Matched: {}", score, matched);
            }
            return (score, 1, matched);
        } else if debug {
            println!("Variant not found in effect weights");
        }
    } else if debug {
        println!("Invalid chromosome or position");
    }
    (0.0, 1, 0)
}
