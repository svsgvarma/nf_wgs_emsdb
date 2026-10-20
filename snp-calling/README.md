# SNP-calling Nextflow pipeline (HPC-portable)

Short-read WGS SNP calling: alignment (BWA / Bowtie2) → samtools sort / add read groups / index →
GATK HaplotypeCaller (GVCF) → GATK joint genotyping → MultiQC.
Cloud/infra code (AWS Batch, S3, Terraform) has been removed; the pipeline runs on any filesystem.

## Layout

```
snp-calling/
├── genomics_snp-calling.nf          # main workflow (BAM list input)
├── genomics_snp-calling_hgtest.nf   # test workflow (FASTQ pairs, GRCh38 GIAB reference)
├── nextflow.config                  # params + executor/container profiles
├── modules/
│   ├── bwa/  bowtie/               # indexing + alignment
│   ├── samtools/index/             # convert, sort, add read groups, index
│   ├── gatk/haplotypecaller/  gatk/jointgenotyping/
│   └── multiqc/
├── Dockerfile                       # local testing only
└── README.md
```

## Input layout

`--input_dir` (default `./input`) must contain:

```
data/sample_bams.txt                 # main workflow: one BAM path per line
data/fastq/<sample>_R{1,2}_*.fastq.gz # hgtest workflow
data/ref/ref.fasta (+ .fai, .dict, intervals.bed)   # main workflow
data/ref/GRCh38_GIAB/...                            # hgtest workflow (names set in the script)
```

For the hgtest workflow the BWA index (`.amb .ann .bwt .pac .sa`) must exist next to the reference
(`bwa index ref.fasta`). Reference paths are defined at the top of each `.nf` file; edit them for your data.

Outputs go to `--outdir` (default `./results`), including `results/reports/` (report, trace, timeline, DAG).

## Run on an HPC

Use the scheduler profile that matches your site, combined with either containers or modules.

```bash
# Slurm + Singularity/Apptainer
nextflow run genomics_snp-calling.nf -profile slurm,singularity --input_dir /scratch/me/project

# Slurm + tools from environment modules (no containers)
module load nextflow bwa samtools gatk multiqc   # names/versions vary per site
nextflow run genomics_snp-calling.nf -profile slurm,nocontainer --input_dir /scratch/me/project
```

Profiles in [nextflow.config](nextflow.config):

| Type | Profiles |
|------|----------|
| Executor | `local`, `slurm`, `pbs`, `sge`, `lsf` |
| Containers | `singularity`, `apptainer`, `docker`, `nocontainer` |

Site adjustments: uncomment/set `process.queue` and `process.clusterOptions` in the chosen profile, set
`NXF_SINGULARITY_CACHEDIR` to a shared location, and adapt the module names to your cluster.
Per-process `cpus`/`memory`/`time` are set in the modules (e.g. 16–64 GB for samtools/GATK steps); tune them
to your nodes. Container images are assigned per process in the `process { withName: ... }` block.

If a module's `publishDir`/work paths should live on scratch, pass `-w /scratch/.../work`.

## Local test with Docker

The [Dockerfile](Dockerfile) bundles Nextflow, Java 17, bwa, bowtie2, samtools, GATK 4.5.0.0 and MultiQC so
the pipeline can be smoke-tested on a workstation. It is not meant for HPC deployment.

```bash
docker build -t rapeseed-snp-calling .

docker run --rm \
  -v "$PWD/input:/pipeline/input" \
  -v "$PWD/results:/pipeline/results" \
  -v "$PWD/work:/pipeline/work" \
  rapeseed-snp-calling genomics_snp-calling.nf --input_dir /pipeline/input
```

PowerShell: replace `$PWD` with `${PWD}`. The container runs with `-profile local,nocontainer`
(tools come from the image). Use a small subset of data; the WGS memory settings (up to 64 GB) may need
lowering for a laptop, e.g. via `-c local_overrides.config`.

## Notes

- `genomics_snp-calling.nf` ends with a commented-out `ALIGN` process template; BAM inputs are used directly.
- Some comments in the scripts/modules still mention S3/Batch from the original cloud setup; they do not
  affect execution.
