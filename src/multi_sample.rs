use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom, BufRead, BufReader};
use flate2::read::GzDecoder;
use std::time::Instant;

const BUFFER_SIZE: usize = 1024 * 1024; // 1MB buffer

trait ReadOnly: Read {}
impl<T: Read> ReadOnly for T {}

trait ReadSeek: Read + Seek {}
impl<T: Read + Seek> ReadSeek for T {}

struct GzDecoderWrapper<R: Read>(GzDecoder<R>);

impl<R: Read> Read for GzDecoderWrapper<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.0.read(buf)
    }
}

impl<R: Read> ReadOnly for GzDecoderWrapper<R> {}

struct ResilienceReader<R: ReadOnly> {
    inner: R,
    buffer: Vec<u8>,
    pos: usize,
    cap: usize,
}

impl<R: ReadOnly> ResilienceReader<R> {
    fn new(inner: R) -> Self {
        Self {
            inner,
            buffer: vec![0; BUFFER_SIZE],
            pos: 0,
            cap: 0,
        }
    }

    fn fill_buffer(&mut self) -> io::Result<()> {
        self.pos = 0;
        self.cap = 0;
        while self.cap < self.buffer.len() {
            match self.inner.read(&mut self.buffer[self.cap..]) {
                Ok(0) => break,
                Ok(n) => self.cap += n,
                Err(ref e) if e.kind() == io::ErrorKind::Interrupted => continue,
                Err(e) => return Err(e),
            }
        }
        Ok(())
    }

    fn read_until(&mut self, delimiter: u8, buf: &mut Vec<u8>) -> io::Result<usize> {
        let mut read = 0;
        loop {
            if self.pos >= self.cap {
                self.fill_buffer()?;
                if self.cap == 0 {
                    return Ok(read);
                }
            }
            
            let (done, used) = {
                let available = &self.buffer[self.pos..self.cap];
                let mut memchr = available.iter().enumerate();
                
                match memchr.find(|&(_, &b)| b == delimiter) {
                    Some((i, _)) => {
                        buf.extend_from_slice(&available[..=i]);
                        (true, i + 1)
                    }
                    None => {
                        buf.extend_from_slice(available);
                        (false, available.len())
                    }
                }
            };
            self.pos += used;
            read += used;
            if done || used == 0 {
                return Ok(read);
            }
        }
    }
}

impl<R: ReadOnly> Read for ResilienceReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.pos >= self.cap {
            self.fill_buffer()?;
            if self.cap == 0 {
                return Ok(0);
            }
        }
        let amt = std::cmp::min(buf.len(), self.cap - self.pos);
        buf[..amt].copy_from_slice(&self.buffer[self.pos..self.pos + amt]);
        self.pos += amt;
        Ok(amt)
    }
}

pub fn calculate_polygenic_score_multi(path: &str, effect_weights: &HashMap<(u8, u32), f32>, debug: bool) -> io::Result<(f64, usize, usize)> {
    let start_time = Instant::now();
    
    if debug {
        println!("Opening file: {}", path);
    }
    let file = File::open(path)?;
    
    let reader: Box<dyn ReadOnly> = if path.ends_with(".gz") {
        if debug {
            println!("Detected gzipped file, using GzDecoder with ResilienceReader");
        }
        Box::new(ResilienceReader::new(GzDecoderWrapper(GzDecoder::new(file))))
    } else {
        if debug {
            println!("Using standard ResilienceReader");
        }
        Box::new(ResilienceReader::new(file))
    };

    if debug {
        println!("Searching for VCF data start...");
    }
    let (header, sample_count) = find_vcf_start(reader)?;

    if debug {
        println!("VCF data start found. Header: {}", header);
        println!("Sample count: {}", sample_count);
        println!("Processing variants...");
    }

    let mut reader = header.reader;
    let mut line_buffer = Vec::new();
    let mut total_score = 0.0;
    let mut total_variants = 0;
    let mut total_matched = 0;

    loop {
        line_buffer.clear();
        match reader.read_until(b'\n', &mut line_buffer) {
            Ok(0) => break, // End of file
            Ok(_) => {
                if !line_buffer.starts_with(b"#") {
                    let line = String::from_utf8_lossy(&line_buffer);
                    let (score, variants, matched) = process_line(&line, effect_weights, sample_count, debug);
                    total_score += score;
                    total_variants += variants;
                    total_matched += matched;

                    if debug && total_variants % 100_000 == 0 {
                        println!("Processed {} variants", total_variants);
                    }
                }
            }
            Err(e) => {
                eprintln!("Error reading line: {}. Skipping to next line.", e);
                // Attempt to recover by seeking to the next newline
                reader.read_until(b'\n', &mut Vec::new())?;
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

    Ok((total_score / sample_count as f64, total_variants, total_matched))
}

struct HeaderInfo {
    header: String,
    sample_count: usize,
    reader: Box<dyn ReadOnly>,
}

fn find_vcf_start(mut reader: Box<dyn ReadOnly>) -> io::Result<HeaderInfo> {
    let mut line = Vec::new();
    let mut last_header_line = String::new();

    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "VCF header not found"));
        }

        let line_str = String::from_utf8_lossy(&line);
        if line_str.starts_with("#CHROM") {
            let sample_count = line_str.split_whitespace().count() - 9;
            return Ok(HeaderInfo {
                header: last_header_line,
                sample_count,
                reader,
            });
        } else if line_str.starts_with('#') {
            last_header_line = line_str.trim().to_string();
        } else {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid VCF format"));
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
            
            if debug && matched > 0 {
                println!("Matched variant at Chr {}, Pos {}. Score: {}, Matched: {}", chr, pos, score, matched);
            }
            
            return (score, 1, matched);
        }
    }
    
    (0.0, 1, 0)
}
