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

    let header = reader.read_header()?;

    for result in reader.records(&header) {
        let record = result?;
        total_variants += 1;

        if let (Some(chr), Some(pos)) = (record.chromosome().to_string().parse::<u8>().ok(), Some(record.position().value() as u32)) {
            if let Some(&weight) = effect_weights.get(&(chr, pos)) {
                if let Some(genotype) = record.genotypes().get(0) {
                    let allele_count = count_alt_alleles(genotype);
                    score += f64::from(weight) * allele_count as f64;
                    matched_variants += 1;
                }
            }
        }
    }

    Ok((score, total_variants, matched_variants))
}

fn count_alt_alleles(genotype: &noodles::vcf::record::genotypes::sample::Value) -> u32 {
    genotype.iter().filter(|&allele| allele == &noodles::vcf::record::genotypes::sample::value::Genotype::Phased(1)).count() as u32
}

pub fn debug_first_lines(path: &str, num_lines: usize) -> io::Result<()> {
    let file = File::open(path)?;
    let buf_reader = BufReader::new(file);
    let mut reader = vcf::Reader::new(buf_reader);

    println!("VCF header:");
    let header = reader.read_header()?;
    println!("{}", header);

    println!("\nFirst {} data lines:", num_lines);
    for (i, result) in reader.records(&header).take(num_lines).enumerate() {
        match result {
            Ok(record) => println!("Line {}: {:?}", i, record),
            Err(e) => println!("Error reading line {}: {:?}", i, e),
        }
    }

    Ok(())
}
