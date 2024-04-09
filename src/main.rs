use std::io::BufReader;
use std::fs::File;
use std::collections::HashMap;

use clap::{Parser};

use rust_htslib::{bam, bam::Read, bam::record::Cigar};
use rust_htslib::bam::ext::BamRecordExtensions;

use bio::io::bed;
use bio::data_structures::interval_tree::*;
use bio::bio_types::genome::AbstractInterval;



#[derive(Parser)]
#[command(name = "seq2c-rs")]
#[command(version)]
#[command(about = "Counts bam coverage of a bed file", long_about = None)]
struct Cli {
    #[arg(short='b', long)]
    bam: String,
    #[arg(short='N', long)]
    sample_name: String,
    #[arg(short='p',long)]
    bed: String,
    #[arg(long)]
    mimic_perl_output: bool,
}


#[derive(Debug, Clone)]
struct RegionWithName {
    name: String,
    count: i64,
}

#[derive(Debug, PartialOrd, Ord, PartialEq, Eq)]
struct OutputRegion {
    name: String,
    start: i64,
    end: i64,
    count: i64
}


fn calculate_coverage(a: std::ops::Range<i64>, b: std::ops::Range<i64>) -> i64 {
    // Find the start and end of the intersection
    let intersection_start = std::cmp::max(a.start, b.start);
    let intersection_end = std::cmp::min(a.end, b.end);
    
    intersection_end - intersection_start + 1
}


fn main(){
    let cli = Cli::parse();
    let sample_name = cli.sample_name;
    let mimic_perl_output = cli.mimic_perl_output;

    let mut bed_map: HashMap<String,IntervalTree<i64, RegionWithName>> = HashMap::new();


    let mut bed_chrom_order = Vec::new();
    let mut reader = File::open(cli.bed).map(BufReader::new).map(bed::Reader::new).unwrap();
    for record in reader.records() {
        let rec = record.expect("Error reading record.");
        //println!("{rec:?}");
        let bed_map_chrom = bed_map.entry(rec.chrom().to_string()).or_default();
        bed_map_chrom.insert(rec.start() as i64..rec.end() as i64,RegionWithName{ name:rec.name().expect("BED record does not define name").to_string(), count: 0});
        bed_chrom_order.push(rec.chrom().to_string());
    }

    // Cleanup chrom ordering from duplicates
    bed_chrom_order.dedup();


    let mut bam = bam::Reader::from_path(cli.bam).unwrap();


    for r in bam.records() {
        let record = r.unwrap();
        if record.tid() < 0 {
            continue;
        }
        if record.is_supplementary() { //skip supplementary aligments
            continue;
        }
        let start = record.reference_start() + 1;  //becuase start position will be included
        let chrom = record.contig();
        let end = start - 1
            + record.cigar()
                .iter()
                .filter_map(|a| match a { 
                   Cigar::Match(l) => Some(l),
                   Cigar::Del(l) => Some(l),
                   _ => None,
                })
                .sum::<u32>() as i64;

        let bed_map_chrom = match bed_map.get_mut(chrom) {
            Some(bed_chrom) => bed_chrom,
            _ => continue,
        };
        for mut interval in bed_map_chrom.find_mut(start-1..end+1) {
            let (interval_start, interval_end) = (interval.interval().start, interval.interval().end);

            let interval_data = interval.data(); //Mutable reference here
            if interval_data.name != "." { //Skip calculation of coverage for unnamed regions
                interval_data.count += calculate_coverage(start..end, interval_start..interval_end);
            }

        }
    }

    // Write the header
    println!("Sample\tGene\tChr\tStart\tEnd\tTag\tLength\tMeanDepth");
    for chrom in bed_chrom_order {
        let chrom_tree = bed_map.get_mut(&chrom).unwrap(); //Safe to unwrap since it's guaranteed that we will have a hit
        let mut output = chrom_tree.find(0..i64::MAX)
                                .map(|e| OutputRegion{start:e.interval().start,
                                                        end:e.interval().end,
                                                    name:e.data().name.clone(),
                                                    count:e.data().count})
                                .collect::<Vec<_>>();

        output.sort();

        let mut current_gene = "";
        let mut total_length = 0;
        let mut current_start = i64::MAX;
        let mut current_end = 0;
        let mut total_count = 0i64;
        
        for region in output.iter() {
            if region.name != current_gene {
                if !current_gene.is_empty() {
                    // Calculate and write aggregated data for the previous gene
                    let mean_depth = if total_length > 0 { total_count as f64 / total_length as f64} else { 0.0 };
                    println!("{sample_name}\t{current_gene}\t{chrom}\t{current_start}\t{current_end}\tWhole-Gene\t{total_length}\t{mean_depth:.2}");
                }
                // Reset
                current_gene = &region.name;
                total_length = 0;
                total_count = 0;
                current_start = region.start;
                current_end = 0;
            }

            // Process current region
            let length = if mimic_perl_output {
                region.end - region.start + 1 //Length in perl version of seq2c calculated +1
            } else {
                region.end - region.start
            };

            let count = region.count;
            if region.end > current_end {
                current_end = region.end;
            }
            println!("{sample_name}\t{}\t{chrom}\t{}\t{}\tAmplicon\t{}\t{:.2}", region.name, region.start, region.end, length, count as f64 /length as f64);
            total_length += length;
            total_count += count;

        }

        // at the end of the vector, write aggregated line for the last gene
        let mean_depth = if total_length > 0 { total_count as f64 / total_length as f64 } else { 0.0 };
        println!("{sample_name}\t{current_gene}\t{chrom}\t{current_start}\t{current_end}\tWhole-Gene\t{total_length}\t{mean_depth:.2}");

    }
}
