use clap::Parser;
use std::collections::{HashMap,HashSet};

#[derive(Debug,Clone)]
pub struct Region{
    pub chr: String,
    pub start: u64,
    pub end: Option<u64>
}

/// Tool that, given a list of variants, will select the ones heterozygous in the DNA bam file, and count reference and alternative read counts in the RNA bam file.
#[derive(Parser,Debug)]
pub struct Config{
    /// Path to the DNA bam file.
    #[arg(long)]
    pub dna:String,

    /// Path to the RNA bam file.
    #[arg(long)]
    pub rna: String,

    /// Path to the vcf file containing variants to consider.
    #[arg(long)]
    pub vcf:String,

    /// Path to the output file.
    #[arg(short,long)]
    pub output:String,

    /// Regions to process (if not provided, will process all regions).
    #[arg(short,long, value_parser=parse_region,num_args =1..)]
    pub regions: Option<Vec<Region>>,

    /// Number of threads to use. If not specified, will use all available threads.
    #[arg(short,long)]
    pub threads: Option<usize>,

    /// Minimum MAPQ for an alignment to be considered.
    #[arg(long,default_value_t = 20)]
    pub min_mapq: u8,

    /// Minimum base quality for a base to be considered.
    #[arg(long,default_value_t = 20)]
    pub min_basequal: u8,

    /// Specify this option to keep duplicate reads. 
    #[arg(long,action)]
    pub keep_duplicates: bool,
    
    /// Specify this option to only also consider indels in addition to SNVs. 
    #[arg(long,action)]
    pub indels: bool,

    /// Minimum number of reads supporting the alternative and reference alleles in the DNA bam file for a variant to be considered.
    #[arg(long,default_value_t = 5)]
    pub min_allele_count_dna: u32,

    /// Minimum variant allele frequency of the minor allele in the DNA bam file for a variant to be considered.
    #[arg(long,default_value_t = 0.30)]
    pub min_vaf_dna: f32,

    /// Minimum coverage in the RNA for a variant to be considered.
    #[arg(long,default_value_t = 6)]
    pub min_coverage_rna: u32,

    /// Maximum number of variants to include in a batch (batches are processed in parallel).
    #[arg(long,default_value_t = 20000)]
    pub max_batchsize: u32
}

fn parse_region(s: &str) -> Result<Region,String>{
    let split :Vec<&str>= s.split(":").collect();
    if split.is_empty() {return Err(String::from("No arguments were provided for the region"));}
    let chr=String::from(split[0]);
    if split.len()==1{
        return Ok(Region{chr,start:0,end:None})
    }
    let split_coords: Vec<&str>=split[1].split("-").collect();
    if split_coords.len()<2 {return Err(String::from("Region must be in the form chr or chr:start-end"));}
    let start=split_coords[0].parse::<u64>();
    if start.is_err() {return Err(format!("Failed to parse as integer: {}",split_coords[0]));}
    let start=start.unwrap();
    let end=split_coords[1].parse::<u64>();
    if end.is_err() {return Err(format!("Failed to parse as integer: {}",split_coords[1]));}
    let end=end.unwrap();
    Ok(Region{chr,start,end:Some(end)})

}


#[derive(Debug,PartialEq)]
pub enum VariantType{
    Snv,
    Insertion,
    Deletion
}

#[derive(Debug)]
pub struct Variant{
    pub chr: String,
    pub pos: u32,
    pub id: String,
    pub variant_type: VariantType,
    pub reference: Vec<u8>,
    pub alternative: Vec<u8>, 
}

impl PartialEq for Variant{
    fn eq(&self, other: &Self) -> bool {
        self.chr == other.chr && self.pos==other.pos && self.reference==other.reference && self.alternative==other.alternative
    }
}



#[derive(Debug)]
pub struct VariantCounter{
    pub variant: Variant,
    pub ref_count_dna: u32,
    pub alt_count_dna: u32,
    pub ref_count_rna: u32,
    pub alt_count_rna: u32,
    pub fragments_processed: HashSet<u64> //qname of reads already counted for this variant, to avoid counting both pairs of an overlapping read pair for one variant.
}

// A list of variants in a region. Batches will be processed in parallel.
// In order to quickly test which variants might overlap a read, store this list of variants as a hashmap of vectors,
// where the key of the hashmap is the position of the variants divided by 10000.
pub struct BatchOfVariants{
    pub chr: String,
    pub start: u32,
    pub end: u32,
    pub hashmap_variantcounters: HashMap<u32,Vec<VariantCounter>>
}