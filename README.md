# fast_ase
Quickly get allelic read counts to estimate allele-specific expression (ASE).

- Takes as input:
  - VCF file of common SNPs (typically restricted to SNPs in genes). Must be indexed.
  - BAM file of whole genome sequencing. Must be indexed.
  - BAM file of RNA-seq. Must be indexed
- Identify the variants in the VCF file that are heterozygous in the DNA.
- For those variants, count the numer of reads supporting the reference and alternative alleles (if the locus is covered).
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
      --dna <DNA>
          Path to the DNA bam file
      --rna <RNA>
          Path to the RNA bam file
      --vcf <VCF>
          Path to the vcf file containing variants to consider
  -o, --output <OUTPUT>
          Path to the output file
  -r, --regions <REGIONS>...
          Regions to process (if not provided, will process all regions)
  -t, --threads <THREADS>
          Number of threads to use. If not specified, will use all available threads
      --min-mapq <MIN_MAPQ>
          Minimum MAPQ for an alignment to be considered [default: 20]
      --min-basequal <MIN_BASEQUAL>
          Minimum base quality for a base to be considered [default: 20]
      --keep-duplicates
          Specify this option to keep duplicate reads
      --indels
          Specify this option to only also consider indels in addition to SNVs
      --min-allele-count-dna <MIN_ALLELE_COUNT_DNA>
          Minimum number of reads supporting the alternative and reference alleles in the DNA bam file for a variant to be considered [default: 5]
      --min-vaf-dna <MIN_VAF_DNA>
          Minimum variant allele frequency of the minor allele in the DNA bam file for a variant to be considered [default: 0.3]
      --min-coverage-rna <MIN_COVERAGE_RNA>
          Minimum coverage in the RNA for a variant to be considered [default: 6]
      --max-batchsize <MAX_BATCHSIZE>
          Maximum number of variants to include in a batch (batches are processed in parallel) [default: 20000]
  -h, --help
          Print help
```

For the VCF file of common SNPs, you can use [this file for hg19](https://drive.google.com/drive/folders/1_Hj7F-13LHz_o8QpU9nOaDvMJdY4n1eZ?usp=drive_link) or [that file for hg38](https://drive.google.com/drive/folders/1-pxEDiml3kQZC7LDbbSnJ0O6BJGbZ4rQ?usp=drive_link). They were downloaded from https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/ and subsetted to SNPs in genes.

## Workflow

fast_ase is integrated into the nextflow workflow [CompEpigen/wf_WGS](https://github.com/CompEpigen/wf_WGS). This can make it easier to run this tool for a large number of samples.
