use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use flate2::read::GzDecoder;


pub fn calculate_polygenic_score_multi(path: &str, effect_weights: &HashMap<(u8, u32), f32>, debug: bool) -> io::Result<(f64, usize, usize)> {
    if debug {
        println!("Opening file: {}", path);
    }
    let file = File::open(path)?;
    let reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        if debug {
            println!("Detected gzipped file, using GzDecoder");
        }
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        if debug {
            println!("Using standard BufReader");
        }
        Box::new(BufReader::new(file))
    };

    let mut lines = reader.lines();
    if debug {
        println!("Searching for header...");
    }
    let header = find_header(&mut lines)?;
    let sample_count = header.split('\t').count() - 9;

    if debug {
        println!("Header found: {}", header);
        println!("Sample count: {}", sample_count);
        println!("Attempting to read first few lines after header:");
    }
    
    let mut first_lines = Vec::new();
    for (i, line) in lines.by_ref().take(5).enumerate() {
        match line {
            Ok(l) => {
                if debug {
                    println!("Line {}: {}", i + 1, l);
                }
                first_lines.push(l);
            },
            Err(e) => {
                if debug {
                    println!("Error reading line {}: {}", i + 1, e);
                }
                return Err(e);
            }
        }
    }
    
    if debug {
        println!("End of first few lines");
        if first_lines.is_empty() {
            println!("Warning: No lines read after header!");
        }
    }

    let mut line_count = 0;
    let mut total_score = 0.0;
    let mut total_variants = 0;
    let mut total_matched = 0;

    // Process the first few lines we've already read
    for line in first_lines {
        let (score, variants, matched) = process_line(&line, effect_weights, sample_count, debug);
        total_score += score;
        total_variants += variants;
        total_matched += matched;
        line_count += 1;
    }

    // Process the rest of the lines
    for line in lines {
        match line {
            Ok(line) => {
                if debug && line_count % 100_000 == 0 {
                    println!("Processing line {}", line_count);
                }
                let (score, variants, matched) = process_line(&line, effect_weights, sample_count, debug);
                total_score += score;
                total_variants += variants;
                total_matched += matched;
                line_count += 1;
            },
            Err(e) => {
                if debug {
                    println!("Error reading line: {}", e);
                }
                // Decide whether to break or continue based on the error
                break;
            }
        }
    }

    if debug {
        println!("Total lines processed: {}", line_count);
    }

    Ok((total_score / sample_count as f64, total_variants, total_matched))
}


fn find_header<B: BufRead>(lines: &mut std::io::Lines<B>) -> io::Result<String> {
    for line in lines.by_ref() {
        let line = line?;
        if line.starts_with("#CHROM") {
            return Ok(line);
        }
    }
    Err(io::Error::new(io::ErrorKind::InvalidData, "Header not found"))
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
            return (score, 1, matched);
        }
    }
    (0.0, 1, 0)
}
