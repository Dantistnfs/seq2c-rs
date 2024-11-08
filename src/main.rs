use std::io::BufReader;
use std::fs::File;
use std::thread::available_parallelism;

use clap::{Parser};

use rust_htslib::{bam, bam::Read, bam::record::Cigar};
use rust_htslib::bam::ext::BamRecordExtensions;

use bio::io::bed;
use bio::bio_types::genome::AbstractInterval;

use std::cell::RefCell;

use coitrees::*;
use rustc_hash::FxHashMap;
use indexmap::IndexMap;
use fnv::FnvBuildHasher;

type FnvIndexMap<K, V> = IndexMap<K, V, FnvBuildHasher>;


#[derive(Parser)]
#[command(name = "seq2c-rs")]
#[command(version)]
#[command(about = "Counts bam coverage of a bed file", long_about = None)]
struct Cli {
    #[arg(short='b', long, help="path to the bam file")]
    bam: String,
    #[arg(short='N', long, help="file name to use in output file")]
    sample_name: String,
    #[arg(short='p',long, help="path to the bed file")]
    bed: String,
    #[arg(long, default_value="true", help="(default: true) enable outputting fragment length - 1, same as perl version of seq2c")]
    mimic_perl_output: bool,
    #[arg(long="threads",default_value="0",help="number of threads to use for bam/cram decompression, default 0 = automatically detect number of cores")]
    threads: usize,
}


#[derive(Debug, Clone)]
struct RegionWithName {
    name: String,
    count: RefCell<i64>,
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



fn update_node(start: i64, end: i64, interval: &IntervalNode<RegionWithName, u32>) {
    let metadata = &interval.metadata;
    if metadata.name != "." { //Skip calculation of coverage for unnamed regions
        let mut count = metadata.count.borrow_mut(); //Mutable borrow, but happens only in one thread, so it's fine
        *count += calculate_coverage(start..end, interval.first as i64..interval.last as i64);
    }
}



fn main(){
    let cli = Cli::parse();
    let sample_name = cli.sample_name;
    let mimic_perl_output = cli.mimic_perl_output;
    eprintln!("Started");

    let bam_threads = if cli.threads == 0 {
            available_parallelism().expect("Wasn't able to automatically reconize number of threads, please set it by setting --threads argument manually").get()
        } else {
            cli.threads
    };
    eprintln!("Using {bam_threads} threads for reading bam file");

    let mut nodes: FxHashMap<String, Vec<Interval<RegionWithName>>> = FxHashMap::default();
    let mut bed_map: FxHashMap<String, COITree<RegionWithName, u32>> = FxHashMap::default();

    eprintln!("Reading bed file");
    let mut bed_chrom_order = Vec::new();
    let mut reader = File::open(cli.bed).map(BufReader::new).map(bed::Reader::new).unwrap();
    for record in reader.records() {
        let rec = record.expect("Error reading record.");
        let node_vec = nodes.entry(rec.chrom().to_string()).or_default();
        node_vec.push(
                        Interval::new(rec.start() as i32, 
                                        rec.end() as i32,
                                        RegionWithName{ 
                                            name:rec.name().expect("BED record does not define name").to_string(), 
                                            count: RefCell::new(0)
                                        }
                                    )
                        );
        bed_chrom_order.push(rec.chrom().to_string());
    }

    for (chrom, chrom_nodes) in nodes {
        bed_map.insert(chrom, COITree::new(&chrom_nodes));
    }
    eprintln!("Reading bed file finished");

    // Cleanup chrom ordering from duplicates
    bed_chrom_order.dedup();

    // Convert COITree to Querent that stores info about last region to optinize serach
    let mut querents = FnvIndexMap::<String, COITreeSortedQuerent<RegionWithName, u32>>::default();
    for (seqname, tree) in &bed_map {
        querents.insert(seqname.clone(), COITreeSortedQuerent::new(tree));
    }

    eprintln!("Starting processing bam file");

    let mut bam = bam::Reader::from_path(cli.bam).unwrap();
    bam.set_threads(bam_threads).expect("Error in setting number of threads for loading bam file");

    for r in bam.rc_records() {
        let record = r.expect("Failure parsing Bam file");
        if record.is_supplementary() { //skip supplementary aligments
            continue;
        }
        if record.tid() < 0 {
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

        let querent_chrom = match querents.get_mut(chrom) {
            Some(querent_chrom) => querent_chrom,
            _ => continue,
        };
        querent_chrom.query((start-1) as i32, (end+1) as i32, |node| {update_node(start, end, node)}); // Runs update_node on
        // each interval in tree that has intersection with query interval
    }

    eprintln!("Finished processing bam file");

    eprintln!("Outputing result into stdout");

    // Prepare the header
    let mut output_string = String::from("Sample\tGene\tChr\tStart\tEnd\tTag\tLength\tMeanDepth\n");

    for chrom in bed_chrom_order {
        let chrom_tree = querents.get_mut(&chrom).unwrap(); //Safe to unwrap since it's guaranteed that we will have a hit
        //let mut output = chrom_tree.iter()
        let mut output = Vec::new();
        chrom_tree.query(0, i32::MAX, |node| {output.push(OutputRegion{start:node.first as i64,
                                                        end:node.last as i64,
                                                        name:node.metadata.name.clone(),
                                                        count:*node.metadata.count.borrow()})
                                            });

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
                    output_string += format!("{sample_name}\t{current_gene}\t{chrom}\t{current_start}\t{current_end}\tWhole-Gene\t{total_length}\t{mean_depth:.2}\n").as_str();
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
            output_string += format!("{sample_name}\t{}\t{chrom}\t{}\t{}\tAmplicon\t{}\t{:.2}\n", region.name, region.start, region.end, length, count as f64 /length as f64).as_str();
            total_length += length;
            total_count += count;

        }

        // at the end of the vector, write aggregated line for the last gene
        let mean_depth = if total_length > 0 { total_count as f64 / total_length as f64 } else { 0.0 };
        output_string += format!("{sample_name}\t{current_gene}\t{chrom}\t{current_start}\t{current_end}\tWhole-Gene\t{total_length}\t{mean_depth:.2}\n").as_str();
    }
    print!("{}", output_string);

    eprintln!("Done");
}
