
# Indexing a Reference sequence and annotation

1. First lets make sure we are where we are supposed to be and that the References directory is available.

    ```bash
    cd /share/workshop/alliance_covid/$USER/covid_swift
    mkdir -p slurmout
    ```

1. To align our data we will need the genome (fasta) and annotation (gtf) for Sars-Cov2. There are many places to find them, but we are going to get them from the [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/1798174254).

    We need to first get the url from the ftp for the genome and annotation gtf.

    We will need:

    *   [Genome sequence](ftp://ftp.ensemblgenomes.org/pub/viruses/fasta/sars_cov_2/dna/Sars_cov_2.ASM985889v3.dna.toplevel.fa.gz)

1. We are going to use an aligner called ['Bowtie2'](https://www.nature.com/articles/nmeth.1923) to align the data. Lets take a look at the help docs for star:

    ```bash
    module load bowtie2
    bowtie2 -h
    ```

1. First we need to index the genome for Bowtie2. Lets pull down a script to index the Ensembl version of the Covid genome. We already copied the genome over yesterday when we pulled over the primer sequences.

    ```bash
    cd resources

    FASTA="NC_045512.2.fasta"

    module load bowtie2

    bowtie2-build \
        --threads 8 \
        ${FASTA} \
        ${FASTA}
    ```


    This step will take about 2 minutes.
