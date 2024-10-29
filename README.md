# fast_ase
Quickly get allelic read counts to estimate allele-specific expression (ASE).

- Takes as input:
  - VCF file of common SNPs (typically restricted to SNPs in genes). Must be indexed.
  - BAM file of whole genome sequencing. Must be indexed.
  - BAM file of RNA-seq. Must be indexed
- Identify the variants in the VCF file that are heterozygous in the DNA.
- For those variants, if they are covered in the RNA-seq data, count the numer of reads supporting the reference and alternative alleles.
- Outputs a tsv file with the allelic read counts in DNA and RNA.
  
Output format:

| contig | position | variantID | refAllele | altAllele | refCount | altCount | refCount_DNA | altCount_DNA |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | 
| 1 | 12267999 | rs1061628 | C | T | 207 | 143 | 36 | 36 |

This table can then be used to infer if some genes are expressed from only one allele (e.g. imprinted genes, genes activated by enhancer hijacking...). In particular, fast_ase can be used as a preprocessing step before running [pyjacker](https://github.com/CompEpigen/pyjacker) to detect enhancer hijacking events in a cohort of cancer patients.

fast_ase is similar to running [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) followed by [GATK ASEReadCounter](https://gatk.broadinstitute.org/hc/en-us/articles/360037428291-ASEReadCounter), but 10-20x faster. fast_ase only tests the variants listed in the VCF file provided as input, which saves a lot of time, but misses 5-10% of variants compared to testing all possible variants with the GATK workflow. In addition, fast_ase does not perform local realignment, which saves time but might not always be 100% accurate, especially in case of indels. Consequently, indel detection is disabled by default, but can be enabled with `--indels`, although the read counts might be inaccurate.

## Usage

```
Usage: fast_ase [OPTIONS] --dna <DNA> --rna <RNA> --vcf <VCF> --output <OUTPUT>

Options:
      --dna <DNA>                      Path to the DNA bam file
      --rna <RNA>                      Path to the RNA bam file
      --vcf <VCF>                      Path to the vcf file containing variants to genotype
  -o, --output <OUTPUT>                Path to the output file
  -r, --regions <REGIONS>...           Regions to process (if not provided, will process all regions)
  -t, --threads <THREADS>              Number of threads to use. If not specified, will use all available threads
      --min-mapq <MIN_MAPQ>            Minimum MAPQ for an alignment to be considered [default: 20]
      --min-basequal <MIN_BASEQUAL>    Minimum base quality for a base to be considered [default: 20]
      --keep-duplicates                Specify this option to keep duplicate reads
      --indels                         Specify this option to only also consider indels in addition to SNVs
      --max-batchsize <MAX_BATCHSIZE>  Minimum MAPQ for an alignment to be considered [default: 20000]
  -h, --help                           Print help
```

## Workflow

fast_ase is integrated into the nextflow workflow [CompEpigen/wf_WGS](https://github.com/CompEpigen/wf_WGS). This can make it easier to run this tool for a large number of samples.
