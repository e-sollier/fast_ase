use std::collections::{HashMap,HashSet};
use rust_htslib::{bcf,errors::Error,errors::Result};

use crate::utils::*;


/// If regions are already provided, simply return them.
/// If not, use one region per chromosome.
pub fn get_regions(config: &Config) ->Vec<Region>{
    let bcf = bcf::IndexedReader::from_path(&config.vcf).expect("Failed to read vcf file.");
    let header=bcf::Read::header(&bcf).clone();

    match &config.regions{
        Some(v) =>v.clone(),
        None =>{
            let mut v:Vec<Region> = Vec::new();
            for i in 0..header.contig_count(){
                let chr=String::from(std::str::from_utf8(header.rid2name(i).unwrap()).unwrap());
                if chr!="Y" && chr!="chrY"{ // at most one chrY so no allele-specific expression
                    v.push(Region{chr,start:0,end:None});
                } 
            }
            v
        }
    }
}

/// Parse a vcf file in a region and returns a vector of batches of variants.
/// Batches contain variants which are close to each other and will be processed in parallel. 
/// A new batch is created when there is a large gap between two variants, or when a batch contains too many variants.
pub fn get_batches_of_variants(config: &Config, region: &Region) -> Result<Vec<BatchOfVariants>>{

    let mut list_of_batches=Vec::new();
    let mut batch_chr = String::new();
    let mut batch_start=0;
    let mut batch_end=0;
    let mut hashmap_variantcounters: HashMap<u32,Vec<VariantCounter>> = HashMap::new();
    let mut batch_size=0;

    let mut bcf = bcf::IndexedReader::from_path(&config.vcf).expect("Failed to read vcf file.");
    let header=bcf::Read::header(&bcf).clone();

    let rid=contig2rid(&region.chr,&header)?;
    bcf.fetch(rid,region.start,region.end)?;
    for record_result in bcf::Read::records(&mut bcf) {
        let record = record_result.expect("Failed to read record.");
        for variant in vcf_record2variants(&header,&record){
            if (!config.indels) && variant.variant_type!=VariantType::Snv {continue;}
            let ind = variant.pos/10000;

            // Create a new batch if the current variant is far from the previous one.
            if variant.chr!=batch_chr || batch_end + 10000 < variant.pos || batch_size>=config.max_batchsize{
                if batch_end>0 {
                    list_of_batches.push(BatchOfVariants{chr:batch_chr,start:batch_start,end:batch_end,hashmap_variantcounters});
                }
                batch_chr = variant.chr.clone();
                batch_start=variant.pos;
                hashmap_variantcounters = HashMap::new();
                batch_size=0;
            } 

            batch_end=variant.pos+1;
            hashmap_variantcounters.entry(ind).or_default();

            //Only add the variant if it is new
            let mut variant_is_new=true;
            let variantcounters=hashmap_variantcounters.get_mut(&ind).unwrap();
            for variant2 in variantcounters.iter(){
                if variant==variant2.variant {variant_is_new=false;}
            }
            if variant_is_new{
                variantcounters.push(VariantCounter{variant,ref_count_dna:0,alt_count_dna:0,ref_count_rna:0,alt_count_rna:0,fragments_processed: HashSet::new()});
                batch_size+=1;
            }
        }
        
    }
    //Add last batch
    list_of_batches.push(BatchOfVariants{chr:batch_chr,start:batch_start,end:batch_end,hashmap_variantcounters});
    
    Ok(list_of_batches)
}




fn vcf_record2variants(header: &bcf::header::HeaderView,record: &bcf::record::Record) -> Vec<Variant>{
    let mut variants: Vec<Variant> = Vec::new();
    let alleles = record.alleles();
    if alleles.len()<2 {return variants}
    let rid = record.rid();
    if rid.is_none() {return variants}
    let rid=rid.unwrap();
    let contig_name=header.rid2name(rid).unwrap();
    let contig_name = std::str::from_utf8(contig_name).unwrap();

    for i in 1..alleles.len(){
        if alleles[0].len()==1 && alleles[i].len()==1 {
            //SNV
            let variant=Variant{chr: String::from(contig_name), pos:record.pos() as u32,
                id: String::from_utf8(record.id()).unwrap(),variant_type: VariantType::Snv, reference:alleles[0].to_vec(),alternative:alleles[i].to_vec()};
            variants.push(variant);
    
        }
        else if  alleles[0].len()==1 {
            //Insertion
            let variant=Variant{chr: String::from(contig_name), pos:record.pos() as u32,
                id: String::from_utf8(record.id()).unwrap(),variant_type: VariantType::Insertion, reference:alleles[0].to_vec(),alternative:alleles[i].to_vec()};
            variants.push(variant);
    
        }
        else if  alleles[i].len()==1 {
            //Deletion
            let variant=Variant{chr: String::from(contig_name), pos:record.pos() as u32,
                id: String::from_utf8(record.id()).unwrap(),variant_type: VariantType::Deletion, reference:alleles[0].to_vec(),alternative:alleles[i].to_vec()};
            variants.push(variant);
        }
    }
    variants    
}

fn contig2rid(chr: &String,header: &bcf::header::HeaderView) -> Result<u32>{
    let mut chr_noprefix=chr.as_bytes();
    if chr_noprefix.len()>=3 && &chr_noprefix[0..3]==b"chr" {chr_noprefix=&chr_noprefix[3..];}

    for i in 0..header.contig_count(){
        let mut contig_name=header.rid2name(i)?;
        if contig_name.len()>=3 && &contig_name[0..3]==b"chr" {contig_name=&contig_name[3..];}
        if contig_name==chr_noprefix {return Ok(i);}
    }
    Err(Error::BcfUnknownContig{contig:chr.clone()})
}