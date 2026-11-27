# SNP-calling Nextflow pipeline (HPC-portable)

Short-read WGS SNP calling from paired FASTQ: BWA-MEM → samtools fixmate/sort/markdup → index + alignment QC →
GATK HaplotypeCaller (GVCF, scattered over interval shards) → per-shard GenomicsDBImport + GenotypeGVCFs →
GatherVcfs → cohort QC report. Runs on any filesystem.

## Layout

```
snp-calling/
├── main.nf                          # the workflow
├── nextflow.config                  # params + executor/container profiles
├── bin/                             # make_shards.py, snp_qc_report.py
├── modules/
│   ├── reference/                  # faidx, dict, BWA index; interval sharding
│   ├── bwa/                        # alignment
│   ├── samtools/                   # markdup, index, qc
│   ├── gatk/scatter/               # sharded calling, joint genotyping, gather
│   ├── qc_report/
│   └── bowtie/  multiqc/  gatk/haplotypecaller/  gatk/jointgenotyping/   # not used by main.nf
└── README.md
```

## Inputs

| Parameter | Meaning |
|-----------|---------|
| `--reads_fastq` | Glob for paired FASTQs, e.g. `'/data/fastq/*_R{1,2}_*.fastq.gz'` (quote it). The part before `_R1`/`_R2` becomes the sample name. |
| `--reference` | Uncompressed FASTA. If `.fai`, `.dict` and the BWA index are not next to it, they are built once and cached in `--index_store`. |
| `--intervals` | Optional BED / `.interval_list` / `.list` restricting calling; otherwise all contigs (filter with `--contig_regex`, `--min_contig_length`). |
| `--scatter_count` | Number of interval shards to call in parallel (default 12). |
| `--cohort_name` | Prefix of the joint VCF (`<cohort_name>.joint.vcf.gz`). |

The defaults for `--reads_fastq` and `--reference` point at a GRCh38/GIAB test layout under `--input_dir`;
set both explicitly for real data.

Outputs go to `--outdir` (default `./results`), including `results/reports/` (report, trace, timeline, DAG).

## Run on an HPC

Use the scheduler profile that matches your site, combined with either containers or modules.

```bash
# Slurm + Singularity/Apptainer
nextflow run main.nf -profile slurm,singularity \
    --reference /scratch/me/ref/ref.fasta --reads_fastq '/scratch/me/fastq/*_R{1,2}_*.fastq.gz'

# Slurm + tools from environment modules (no containers)
module load nextflow bwa samtools gatk python   # names/versions vary per site
nextflow run main.nf -profile slurm,nocontainer \
    --reference /scratch/me/ref/ref.fasta --reads_fastq '/scratch/me/fastq/*_R{1,2}_*.fastq.gz'
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

The [Dockerfile](../Dockerfile) in the repo root bundles Nextflow, Java 17, bwa, bowtie2, samtools, GATK 4.5.0.0 and MultiQC so
the pipeline can be smoke-tested on a workstation. It is not meant for HPC deployment.

```bash
# run from the repo root
docker build -t nf-wgs-emsdb .

docker run --rm \
  -v "$PWD/input:/pipeline/input" \
  -v "$PWD/results:/pipeline/snp-calling/results" \
  -v "$PWD/work:/pipeline/snp-calling/work" \
  nf-wgs-emsdb main.nf \
    --reference /pipeline/input/ref/ref.fasta \
    --reads_fastq '/pipeline/input/fastq/*_R{1,2}_*.fastq.gz'
```

PowerShell: replace `$PWD` with `${PWD}`. The container runs with `-profile local,nocontainer`
(tools come from the image). Use a small subset of data; the WGS memory settings (up to 64 GB) may need
lowering for a laptop, e.g. via `-c local_overrides.config`.

## Notes

- Some comments in the scripts/modules still mention S3/Batch from the original cloud setup; they do not
  affect execution.
