# Rapeseed WGS – EMS mutant analysis

Nextflow pipelines for whole-genome sequencing (WGS) analysis of **rapeseed (*Brassica napus*) EMS-mutagenized
lines**. The main goal is to call SNPs/small variants across a cohort of mutant lines, from which EMS-induced
mutations (mainly G/C→A/T transitions) can be identified and filtered downstream.

The code is portable to any HPC (Slurm, PBS, SGE, LSF) or a local machine; no cloud/infra code is included.

## Repository layout

```
rapeseed_wgs/
├── README.md                 # this file
└── snp-calling/              # SNP/variant-calling pipeline (Nextflow DSL2)
    ├── genomics_snp-calling.nf
    ├── genomics_snp-calling_hgtest.nf
    ├── nextflow.config
    ├── modules/
    ├── Dockerfile
    └── README.md             # detailed usage, profiles, input layout
```

## Pipelines

| Script | Starts from | What it does |
|--------|-------------|--------------|
| [genomics_snp-calling.nf](snp-calling/genomics_snp-calling.nf) | Aligned BAMs (list in `sample_bams.txt`) | Index BAMs (samtools) → per-sample GVCF (GATK HaplotypeCaller) → joint genotyping across all samples (GATK GenomicsDBImport + GenotypeGVCFs) → cohort VCF. Main workflow for already-aligned rapeseed lines. |
| [genomics_snp-calling_hgtest.nf](snp-calling/genomics_snp-calling_hgtest.nf) | Paired-end FASTQ | Full workflow: BWA-MEM alignment → samtools sort → index → HaplotypeCaller → joint genotyping. Currently configured with a human GRCh38 / NIST GIAB test dataset to validate the pipeline; point the reference and read paths to the *B. napus* genome and EMS-line FASTQs for real analysis. |

### Modules ([snp-calling/modules](snp-calling/modules))

- `bwa/` – reference indexing and BWA-MEM alignment (streams directly to BAM).
- `bowtie/` – Bowtie2 indexing and alignment (alternative aligner, not wired into the workflows).
- `samtools/index/` – SAM→BAM conversion, sorting, read-group addition, BAM indexing.
- `gatk/haplotypecaller/` – per-sample variant calling in GVCF mode.
- `gatk/jointgenotyping/` – combine GVCFs and jointly genotype the cohort.
- `multiqc/` – QC report aggregation (currently commented out in the workflows).

## Quick start

```bash
cd snp-calling
nextflow run genomics_snp-calling.nf -profile slurm,singularity --input_dir /path/to/project
```

See [snp-calling/README.md](snp-calling/README.md) for the input layout, executor/container profiles, HPC module
usage and the local Docker test.

## Notes for EMS data

- Joint genotyping across all lines (and the wild-type/parent) lets you separate EMS-induced mutations from
  background variation; filtering (e.g. mutation-type, allele fraction, depth, shared-variant removal) is done
  downstream of the cohort VCF and is not part of these pipelines yet.
- *B. napus* is allotetraploid (AACC) with a large, repetitive genome; adjust memory/time in the modules and
  the intervals file to your reference.
