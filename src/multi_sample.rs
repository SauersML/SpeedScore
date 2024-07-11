use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, BufRead, BufReader, Seek, SeekFrom};
use std::time::Instant;
use flate2::read::MultiGzDecoder;
use flate2::Decompress;

const BGZF_MAGIC: &[u8] = &[0x1f, 0x8b, 0x08, 0x04];
const GZIP_MAGIC: &[u8] = &[0x1f, 0x8b, 0x08];

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

#[derive(Debug, PartialEq)]
enum FileFormat {
    Gzip,
    Bgzip,
    Uncompressed,
}

struct BGZFBlock {
    uncompressed_data: Vec<u8>,
}

struct BGZFReader<R: Read + Seek> {
    inner: R,
    current_block: Option<BGZFBlock>,
    buffer_pos: usize,
}

impl<R: Read + Seek> BGZFReader<R> {
    fn new(inner: R) -> Self {
        BGZFReader {
            inner,
            current_block: None,
            buffer_pos: 0,
        }
    }

    fn read_block(&mut self) -> io::Result<bool> {
        let mut header = [0u8; 18];
        if self.inner.read_exact(&mut header).is_err() {
            return Ok(false);
        }

        if &header[0..4] != BGZF_MAGIC {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid BGZF block"));
        }

        let block_size = u16::from_le_bytes([header[16], header[17]]) as usize + 1;
        let mut compressed_data = vec![0u8; block_size - 18];
        self.inner.read_exact(&mut compressed_data)?;

        let mut decompressor = Decompress::new(true);
        let mut uncompressed_data = Vec::with_capacity(65536);
        
        decompressor.decompress_vec(
            &compressed_data,
            &mut uncompressed_data,
            flate2::FlushDecompress::Finish,
        )?;

        self.current_block = Some(BGZFBlock {
            uncompressed_data,
        });
        self.buffer_pos = 0;

        Ok(true)
    }
}

impl<R: Read + Seek> Read for BGZFReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut total_read = 0;

        while total_read < buf.len() {
            if self.current_block.is_none() || self.buffer_pos >= self.current_block.as_ref().unwrap().uncompressed_data.len() {
                if !self.read_block()? {
                    break;
                }
            }

            let block = self.current_block.as_ref().unwrap();
            let remaining = block.uncompressed_data.len() - self.buffer_pos;
            let to_read = std::cmp::min(remaining, buf.len() - total_read);

            buf[total_read..total_read + to_read].copy_from_slice(&block.uncompressed_data[self.buffer_pos..self.buffer_pos + to_read]);
            self.buffer_pos += to_read;
            total_read += to_read;
        }

        Ok(total_read)
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
            self.reader.read_line(&mut line)?;
            if line.starts_with("#CHROM") {
                self.sample_count = line.split_whitespace().count() - 9;
                return Ok(());
            } else if !line.starts_with('#') {
                return Err(VcfError::InvalidFormat("VCF header not found".to_string()));
            }
            line.clear();
        }
    }

    fn read_line(&mut self, buf: &mut String) -> Result<usize, VcfError> {
        buf.clear();
        Ok(self.reader.read_line(buf)?)
    }
}

fn detect_file_format<R: Read + Seek>(mut reader: R) -> io::Result<FileFormat> {
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;
    reader.seek(SeekFrom::Start(0))?;

    if magic.starts_with(BGZF_MAGIC) {
        Ok(FileFormat::Bgzip)
    } else if magic.starts_with(GZIP_MAGIC) {
        Ok(FileFormat::Gzip)
    } else {
        Ok(FileFormat::Uncompressed)
    }
}

fn open_vcf_reader(path: &str) -> Result<Box<dyn Read>, VcfError> {
    let file = File::open(path)?;
    let format = detect_file_format(&file)?;

    match format {
        FileFormat::Bgzip => Ok(Box::new(BGZFReader::new(file))),
        FileFormat::Gzip => Ok(Box::new(MultiGzDecoder::new(file))),
        FileFormat::Uncompressed => Ok(Box::new(file)),
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
        println!("Effect weights loaded: {:?}", effect_weights.iter().take(5).collect::<Vec<_>>());
    }

    let reader = open_vcf_reader(path)?;
    let mut vcf_reader = VcfReader::new(reader)?;

    if debug {
        println!("VCF data start found.");
        println!("Sample count: {}", vcf_reader.sample_count);
        println!("Processing variants...");
    }

    let mut line = String::new();
    let mut total_score = 0.0;
    let mut total_variants = 0;
    let mut total_matched = 0;
    let mut in_data_section = false;

    if debug {
        println!("Starting to read lines...");
    }

    while vcf_reader.read_line(&mut line)? > 0 {
        if debug {
            println!("Line read: {}", line);
        }

        if line.starts_with("#CHROM") {
            in_data_section = true;
            if debug {
                println!("Found #CHROM line: {}", line);
            }
            continue;
        }

        if in_data_section {
            if debug {
                println!("Processing variant line: {}", line);
            }
            let (score, variants, matched) = process_line(&line, effect_weights, vcf_reader.sample_count, debug);
            total_score += score;
            total_variants += variants;
            total_matched += matched;

            if debug && total_variants % 100_000 == 0 {
                println!("Processed {} variants", total_variants);
            }
        } else {
            if debug {
                println!("Skipping header/comment line: {}", line);
            }
        }
    }

    let duration = start_time.elapsed();

    if debug {
        println!("Finished reading lines.");
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
