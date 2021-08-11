
# Indexing a Reference sequence and annotation

1. First lets make sure we are where we are supposed to be and that the References directory is available.

    ```bash
    cd /share/workshop/alliance_covid/$USER/covid_swift
    mkdir -p slurmout
    ```

1. To align our data we will need the genome (fasta) and annotation (gtf) for Sars-Cov2. There are many places to find them, but we are going to get them from the [Ensembl](https://covid-19.ensembl.org/info/data/export.html).

    We need to first get the url from the ftp for the genome and annotation gtf.

    We will need:

    *   [Genome sequence](ftp://ftp.ensemblgenomes.org/pub/viruses/fasta/sars_cov_2/dna/Sars_cov_2.ASM985889v3.dna.toplevel.fa.gz)
    *   [Basic gene annotation (GTF)](ftp://ftp.ensemblgenomes.org/pub/viruses/gtf/sars_cov_2/Sars_cov_2.ASM985889v3.101.gtf.gz)

1. We are going to use an aligner called ['Bowtie2'](https://www.nature.com/articles/nmeth.1923) to align the data. Lets take a look at the help docs for star:

    ```bash
    module load bowtie2
    bowtie2 -h
    ```

1. First we need to index the genome for Bowtie2. Lets pull down a script to index the Ensembl version of the Covid genome.

    ```bash
    wget https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2021-Alliance-Makerere_Covid/master/software_scripts/scripts/bowtie2_index.sh
    less bowtie2_index.sh
    ```

    <div class="script">#!/bin/bash

    start=`date +%s`
    echo $HOSTNAME

    outpath="resources"
    mkdir -p ${outpath}

    cd ${outpath}

    wget ftp://ftp.ensemblgenomes.org/pub/viruses/fasta/sars_cov_2/dna/Sars_cov_2.ASM985889v3.dna.toplevel.fa.gz
    gunzip Sars_cov_2.ASM985889v3.dna.toplevel.fa.gz
    FASTA="Sars_cov_2.ASM985889v3.dna.toplevel.fa"

    wget ftp://ftp.ensemblgenomes.org/pub/viruses/gtf/sars_cov_2/Sars_cov_2.ASM985889v3.101.gtf.gz
    gunzip Sars_cov_2.ASM985889v3.101.gtf.gz
    GTF="Sars_cov_2.ASM985889v3.101.gtf"

    module load bowtie2

    call="bowtie2-build \
        --threads 8 \
        ${FASTA} \
        ${FASTA}"

    echo $call
    eval $call

    end=`date +%s`
    runtime=$((end-start))
    echo $runtime
    </div>

    When you are done, type "q" to exit.

    1. The script uses wget to download the fasta and GTF files from Ensembl using the links you found earlier.
    1. Uncompresses them using gunzip.
    1. Run bowtie2-build to generate indexes



1. Run star indexing when ready.

    ```bash
    bash bowtie2_index.sh
    ```

    This step will take about 2 minutes.
