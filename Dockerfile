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
      bcftools \
      htslib \
      snpeff=5.2 \
      gatk4=4.5.0.0 \
      multiqc \
 && micromamba clean -a -y

RUN SNPEFF_DIR="$(find /opt/conda/share -maxdepth 1 -type d -name 'snpeff-*' | head -n 1)" \
 && test -n "$SNPEFF_DIR" \
 && ln -s "$SNPEFF_DIR" /opt/conda/share/snpeff

ENV PATH="/opt/conda/bin:${PATH}"
ENV NXF_HOME=/tmp/.nextflow
ENV SNPEFF_HOME=/opt/conda/share/snpeff
ENV SNPEFF_JAR=/opt/conda/share/snpeff/snpEff.jar
ENV SNPSIFT_JAR=/opt/conda/share/snpeff/SnpSift.jar

WORKDIR /pipeline
COPY --chown=$MAMBA_USER:$MAMBA_USER . /pipeline
WORKDIR /pipeline/snp-calling

# Pass pipeline script and params at run time, e.g.
#   docker run --rm -v ${PWD}/input:/pipeline/input -v ${PWD}/results:/pipeline/results \
#     nf-wgs-emsdb genomics_snp-calling.nf --input_dir /pipeline/input
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "nextflow", "run", "-profile", "local,nocontainer"]
CMD ["genomics_snp-calling.nf"]
