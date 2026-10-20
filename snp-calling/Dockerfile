# Local test image only: bundles Nextflow + every tool the pipeline calls so it can run
# with `-profile local,nocontainer`. On the HPC, skip this image and load the equivalent
# environment modules (or use -profile slurm,singularity).
FROM mambaorg/micromamba:1.5.10

USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps curl ca-certificates \
 && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

RUN micromamba install -y -n base -c conda-forge -c bioconda \
      nextflow \
      openjdk=17 \
      bwa \
      bowtie2 \
      samtools \
      gatk4=4.5.0.0 \
      multiqc \
 && micromamba clean -a -y

ENV PATH="/opt/conda/bin:${PATH}"
ENV NXF_HOME=/tmp/.nextflow

WORKDIR /pipeline
COPY --chown=$MAMBA_USER:$MAMBA_USER . /pipeline

# Pass pipeline script and params at run time, e.g.
#   docker run --rm -v ${PWD}/input:/pipeline/input -v ${PWD}/results:/pipeline/results \
#     rapeseed-snp-calling genomics_snp-calling.nf --input_dir /pipeline/input
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "nextflow", "run", "-profile", "local,nocontainer"]
CMD ["genomics_snp-calling.nf"]
