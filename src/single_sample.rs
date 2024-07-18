use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader};
use noodles::vcf;

pub fn calculate_polygenic_score(path: &str, effect_weights: &HashMap<(u8, u32), f32>) -> io::Result<(f64, usize, usize)> {
    let file = File::open(path)?;
    let buf_reader = BufReader::new(file);
    let mut reader = vcf::Reader::new(buf_reader);

    let mut score = 0.0;
    let mut total_variants = 0;
    let mut matched_variants = 0;

    reader.read_header()?;

    for result in reader.records() {
        let record = result?;
        total_variants += 1;

        if let (Some(chr), Some(pos)) = (record.chromosome().parse::<u8>().ok(), record.position().get() as u32) {
            if let Some(&weight) = effect_weights.get(&(chr, pos)) {
                if let Some(genotype) = record.genotypes().get(0) {
                    let allele_count = genotype.count();
                    score += f64::from(weight) * allele_count as f64;
                    matched_variants += 1;
                }
            }
        }
    }

    Ok((score, total_variants, matched_variants))
}

pub fn debug_first_lines(path: &str, num_lines: usize) -> io::Result<()> {
    let file = File::open(path)?;
    let buf_reader = BufReader::new(file);
    let mut reader = vcf::Reader::new(buf_reader);

    println!("VCF header:");
    match reader.read_header() {
        Ok(header) => println!("{}", header),
        Err(e) => println!("Error reading header: {:?}", e),
    }

    println!("\nFirst {} data lines:", num_lines);
    for (i, result) in reader.records().take(num_lines).enumerate() {
        match result {
            Ok(record) => println!("Line {}: {:?}", i, record),
            Err(e) => println!("Error reading line {}: {:?}", i, e),
        }
    }

    Ok(())
}
