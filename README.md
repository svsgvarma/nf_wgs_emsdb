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

## Submitting jobs to an HPC cluster

Nextflow runs as a lightweight *head job* that itself submits every pipeline step to the scheduler, so only the
head job needs to be submitted (or run in `tmux`/`screen` on a login node for short tests).

1. **Get the code and Nextflow** on the cluster:
   ```bash
   git clone <repo-url> rapeseed_wgs && cd rapeseed_wgs/snp-calling
   module avail nextflow            # module names/versions differ per site
   module load nextflow java        # Nextflow needs Java 11+ (17 recommended)
   ```
2. **Put data on shared/scratch storage** using the layout in [snp-calling/README.md](snp-calling/README.md)
   (`data/ref/...`, `data/sample_bams.txt` or `data/fastq/...`). Use absolute paths in `sample_bams.txt`.
3. **Set site options** in the profile you use in [snp-calling/nextflow.config](snp-calling/nextflow.config)
   (e.g. `process.queue`, `process.clusterOptions = '--account=<project>'`), or pass a site config with `-c site.config`.
4. **Choose how tools are provided**: `-profile slurm,singularity` (containers) or `-profile slurm,nocontainer`
   after loading the tool modules (`bwa samtools gatk multiqc`) in the job script.
5. **Submit the head job** with a script like `run_nf.sbatch` (Slurm shown; for PBS/SGE/LSF swap the directives
   and use `-profile pbs|sge|lsf`):
   ```bash
   #!/bin/bash
   #SBATCH --job-name=nf-snp-calling
   #SBATCH --partition=<partition>
   #SBATCH --account=<project>
   #SBATCH --cpus-per-task=2
   #SBATCH --mem=8G
   #SBATCH --time=72:00:00
   #SBATCH --output=nf_head_%j.log

   module load nextflow java
   # module load bwa samtools gatk multiqc      # only for -profile ...,nocontainer

   export NXF_OPTS='-Xms512m -Xmx4g'
   export NXF_SINGULARITY_CACHEDIR=/scratch/$USER/singularity_cache   # only for singularity

   nextflow run genomics_snp-calling.nf \
       -profile slurm,singularity \
       --input_dir /scratch/$USER/rapeseed_wgs \
       --outdir    /scratch/$USER/rapeseed_wgs/results \
       -w          /scratch/$USER/rapeseed_wgs/work \
       -resume
   ```
   ```bash
   sbatch run_nf.sbatch
   ```
6. **Monitor and manage**:
   ```bash
   squeue -u $USER                  # head job + submitted task jobs
   tail -f nf_head_<jobid>.log      # pipeline progress
   nextflow log                     # run history; reports are in results/reports/
   scancel <jobid>                  # stop the head job (cancel its task jobs too if needed)
   ```
   Re-submit the same command with `-resume` to continue after a failure or a time-limit kill; completed steps
   are reused from the work directory.

Tips: keep `work/` on fast scratch (not home); request a head-job walltime longer than the whole run; lower or raise
per-process `cpus`/`memory`/`time` in the modules to fit your queue limits; delete `work/` once results are verified.

## Notes for EMS data

- Joint genotyping across all lines (and the wild-type/parent) lets you separate EMS-induced mutations from
  background variation; filtering (e.g. mutation-type, allele fraction, depth, shared-variant removal) is done
  downstream of the cohort VCF and is not part of these pipelines yet.
- *B. napus* is allotetraploid (AACC) with a large, repetitive genome; adjust memory/time in the modules and
  the intervals file to your reference.
