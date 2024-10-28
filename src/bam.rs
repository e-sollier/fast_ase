use rust_htslib::{bam, bam::Read, bam::record::Cigar, errors::Result,errors::Error};
use std::collections::HashMap;
use twox_hash::XxHash64;
use std::hash::{Hash, Hasher};

use crate::utils::*;

/// Given a batch of variants, select the variants heterozygous in the DNA bam file, and count the alleles in the RNA bam file for those heterozygous variants. 
pub fn process_batch(batch: &mut BatchOfVariants, config: &Config) -> Vec<VariantCounter>{

    // Count allelic reads in DNA
    let mut bam = bam::IndexedReader::from_path(&config.dna).expect("Failed to open bam file.");
    let tid = contig2tid(&batch.chr,bam.header()).unwrap();
    let fetch_res=bam.fetch((tid,batch.start as i64,batch.end as i64));
    if fetch_res.is_err() {panic!("Failed to fetch region {}:{}-{} in file {}",&batch.chr,batch.start,batch.end,&config.dna);}
    for r in bam.records() {
        let record = r.unwrap();
        if record.mapq()<config.min_mapq {continue;}
        update_variantcounters_record(&record,&mut batch.hashmap_variantcounters,config,false);
    }

    // Filter to only keep heterozygous variants
    for (_,list) in batch.hashmap_variantcounters.iter_mut(){
        list.retain(|variant_counter| {
            if variant_counter.ref_count_dna>5 && variant_counter.alt_count_dna>5{
                let vaf = (variant_counter.alt_count_dna as f32) / ((variant_counter.alt_count_dna as f32) + (variant_counter.ref_count_dna as f32));
                return (0.30..=0.70).contains(&vaf);
            }
            false
        });
        list.iter_mut().for_each(|variant_counter| variant_counter.fragments_processed.clear());
    }
    batch.hashmap_variantcounters.retain(|_,l| !l.is_empty());

    // Count allelic reads in RNA
    let mut bam_rna = bam::IndexedReader::from_path(&config.rna).unwrap();
    let tid_rna = contig2tid(&batch.chr,bam_rna.header()).unwrap();
    let fetch_res=bam_rna.fetch((tid_rna,batch.start as i64,batch.end as i64));
    if fetch_res.is_err() {panic!("Failed to fetch region {}:{}-{} in file {}",&batch.chr,batch.start,batch.end,&config.rna);}
    for r in bam_rna.records() {
        let record = r.unwrap();
        update_variantcounters_record(&record,&mut batch.hashmap_variantcounters,config,true);
    }

    // Filter to only keep variants with coverage
    for (_,list) in batch.hashmap_variantcounters.iter_mut(){
        list.retain(|variant_counter| {variant_counter.ref_count_rna+variant_counter.alt_count_rna>5});
        list.iter_mut().for_each(|variant_counter| variant_counter.fragments_processed.clear());
    }
    batch.hashmap_variantcounters.retain(|_,l| !l.is_empty());

    // Concatenate all results
    let mut sorted_keys: Vec<u32>= batch.hashmap_variantcounters.keys().cloned().collect();
    sorted_keys.sort();

    let mut res:Vec<VariantCounter> = Vec::new();
    for k in sorted_keys.iter(){
        res.append(batch.hashmap_variantcounters.get_mut(k).unwrap());
    }
    res
}

/// For a bam record, find which variants are covered, and update the allelic read counts for those.
fn update_variantcounters_record(record: &bam::record::Record, hashmap: &mut HashMap<u32,Vec<VariantCounter>>,config: &Config,is_rna: bool){

    // Compute a hash for the record qname, in order to only count one alignment per read pair for each variant.
    let mut hasher = XxHash64::with_seed(0);
    record.qname().hash(&mut hasher);
    let fragment_name_hash = hasher.finish();

    let index_start=(record.pos() as u32) / 10000;
    let index_end = record_endpos(record) / 10000;
    for index in index_start..(index_end+1) {
        let b=hashmap.get_mut(&index);
        if b.is_none() {continue}
        let b = b.unwrap();
        for variant_counter in b.iter_mut(){
            if variant_counter.fragments_processed.contains(&fragment_name_hash) {continue;} // Only count one alignment per read pair, for each variant.
            if let Some(variant_present) = variant_in_record(&variant_counter.variant, record, config){
                if variant_present{
                    if is_rna {variant_counter.alt_count_rna+=1;}
                    else {variant_counter.alt_count_dna+=1;}
                }
                else{
                    if is_rna {variant_counter.ref_count_rna+=1;}
                    else {variant_counter.ref_count_dna+=1;}
                }
                variant_counter.fragments_processed.insert(fragment_name_hash);
            }
        }
    }
}

/// Return Some(true) if the alternative allele is present, Some(false) if the reference allele is present, or None if the variant does not overlap the record.
fn variant_in_record(variant: &Variant, record: &bam::record::Record, config: &Config) -> Option<bool>{
    if record.is_duplicate() && (!config.keep_duplicates) {return None;}
    if record.mapq() < config.min_mapq {return None;}
    match variant.variant_type{
        VariantType::Snv =>{
            test_snv(variant,record,config.min_basequal)
        },
        VariantType::Deletion =>{
            test_deletion(variant,record)
        },
        VariantType::Insertion =>{
            test_insertion(variant,record)
        }
    }
}

/// Test if a particular SNV is present in the record.
fn test_snv(variant: &Variant, record: &bam::record::Record, min_basequal:u8) -> Option<bool>{
    let cigar = &record.cigar();
    let refpos=variant.pos;
    let mut pos=cigar.pos() as u32;
    let mut seqpos=0;
    if pos>refpos {return None;}

    for x in cigar.iter() {
        match x {
            Cigar::Match(length) | Cigar::Equal(length) | Cigar::Diff(length) => {
                if pos + length > refpos {
                    // Found the position in the sequence corresponding to the SNV
                    seqpos += refpos - pos;
                    if record.qual()[seqpos as usize] < min_basequal {
                        return None;
                    }
                    let base = record.seq()[seqpos as usize];
                    if base == variant.reference[0] {
                        return Some(false);
                    } else if base == variant.alternative[0] {
                        return Some(true);
                    } else {
                        return None;
                    }
                }
                pos += length;
                seqpos += length;
            },
            Cigar::Del(length) | Cigar::RefSkip(length) => {
                pos += length;
                if pos > refpos {
                    return None;
                }
            },
            Cigar::Ins(length) | Cigar::SoftClip(length) => {
                seqpos += length;
            },
            Cigar::HardClip(_) | Cigar::Pad(_) => (),
        }
    }
    None
}

/// Test if a particular deletion is present in the record.
fn test_deletion(variant: &Variant, record: &bam::record::Record) -> Option<bool>{
    // Variant starts at pos, but the deletion starts at pos+1.
    let cigar = &record.cigar();
    let mut pos=cigar.pos() as u32;
    let mut seq_pos=0;
    if pos>variant.pos+1 {return None;}

    for x in cigar.iter(){
        match x{
            Cigar::Match(length) | Cigar::Equal(length) | Cigar::Diff(length) =>{
                if pos+length>variant.pos+1 {return Some(false);}
                pos+=length;
                seq_pos+=length;
            },
            Cigar::Del(length)  => {
                if pos==variant.pos+1 && *length==((variant.reference.len()-1) as u32) {return Some(true)}
                pos+=length;
                if pos>variant.pos {return None;}
            },
            Cigar::RefSkip(length)=>{
                pos+=length;
                if pos>variant.pos {return None;}
            }
            Cigar::Ins(length) | Cigar::SoftClip(length)=> {seq_pos+=length;},
            Cigar::HardClip(_) | Cigar::Pad(_) => ()
        }
    }
    None
}

/// Test if a particular insertion is present in the record.
fn test_insertion(variant: &Variant, record: &bam::record::Record) -> Option<bool>{
    // Variant starts at pos, but the insertion starts at pos+1.
    let cigar = record.cigar();
    let mut pos=cigar.pos() as u32;
    let mut seq_pos=0;
    if pos>variant.pos+1 {return None;}

    for x in cigar.iter(){
        match x{
            Cigar::Match(length) | Cigar::Equal(length) | Cigar::Diff(length) =>{
                if pos+length>variant.pos+1 {return Some(false);}
                pos+=length;
                seq_pos+=length;
            },
            Cigar::Del(length) | Cigar::RefSkip(length) => {
                pos+=length;
                if pos>variant.pos {return None;}
            },
            Cigar::Ins(length) => {
                if pos==variant.pos+1 {
                    if seq_pos+length<(record.seq_len() as u32)&& *length ==((variant.alternative.len()-1) as u32){
                        let mut same_sequence=true;
                        for i in 0..*length{
                            if record.seq()[(seq_pos+i) as usize] != variant.alternative[(1+i) as usize]{
                                same_sequence=false;
                                break;
                            }
                        }
                        return Some(same_sequence);
                    }
                    else {return Some(false)}
                }
                seq_pos+=length;
            },
            Cigar::SoftClip(length) => {seq_pos+=length;}
            Cigar::HardClip(_) | Cigar::Pad(_) => ()
        }
    }
    None
}

/// Compute the coordinate in the reference for the rightmost position in the alignment.
fn record_endpos(record: &bam::record::Record) -> u32{
    let cigar = record.cigar();
    let mut pos=cigar.pos() as u32;
    for x in cigar.iter(){
        match x{
            Cigar::Match(length) | Cigar::Equal(length) | Cigar::Diff(length) | Cigar::Del(length) | Cigar::RefSkip(length)=>{
                pos+=length;
            },
            _ => ()
        }
    }
    pos
}

fn contig2tid(chr: &String, header: &bam::HeaderView) -> Result<u32> {
    let mut chr_noprefix=chr.as_bytes();
    if chr_noprefix.len()>=3 && &chr_noprefix[0..3]==b"chr" {chr_noprefix=&chr_noprefix[3..];}
    for i in 0..header.target_count(){
        let mut contig_name=header.tid2name(i);
        if contig_name.len()>=3 && &contig_name[0..3]==b"chr" {contig_name=&contig_name[3..];}
        if contig_name==chr_noprefix {return Ok(i);}
    }
    Err(Error::Fetch)
}