use std::vec::Vec;
use std::fs::File;
use std::io::{BufWriter, Write};

use rayon::prelude::*;
use clap::Parser;


mod utils;
mod vcf;
mod bam;
use utils::*;



// The genome is split into regions (chromosomes by default, or user-specified).
// Then, each region will be split into batches of variants 
// (a new batch is created if there is a large distance between two variants [>10kb by default], or if too many variants are already in the batch [>20000 by default]).
// Regions are processed sequentially, to ensure that the memory footprint is not too large, and batches for each region are processed in parallel.

fn main() {
    let config=Config::parse();

    if let Some(t)=config.threads {
        rayon::ThreadPoolBuilder::new().num_threads(t).build_global().unwrap();
    }

    let regions: Vec<Region> = vcf::get_regions(&config);



    // Output file
    let f = File::create(&config.output).expect("unable to create file");
    let mut f = BufWriter::new(f);
    writeln!(f,"contig\tposition\tvariantID\trefAllele\taltAllele\trefCount\taltCount\trefCount_DNA\taltCount_DNA").expect("unable to write");


    // Parse regions sequentially
    for region in &regions{

        // Parse vcf to find variants in the region
        let region_batches = vcf::get_batches_of_variants(&config,region);
        if let Err(e)=region_batches {
            println!("Failed to process region {:?} in vcf file. {}",&region,e);
            continue;
        }
        let mut region_batches=region_batches.unwrap();
        
        // Process all batches of variants for a region in parallel
        let lists: Vec<Vec<VariantCounter>>= region_batches.par_iter_mut().map(|batch: &mut BatchOfVariants| bam::process_batch(batch, &config)).filter(|v| !v.is_empty()).collect();
        
        // Write to the output file.
        for v in lists{
            for variant_counter in v{
                writeln!(f, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", variant_counter.variant.chr,variant_counter.variant.pos+1,variant_counter.variant.id,String::from_utf8(variant_counter.variant.reference).unwrap(),String::from_utf8(variant_counter.variant.alternative).unwrap(),
            variant_counter.ref_count_rna,variant_counter.alt_count_rna,variant_counter.ref_count_dna,variant_counter.alt_count_dna).expect("unable to write");
            }
        }
    }
    
}
