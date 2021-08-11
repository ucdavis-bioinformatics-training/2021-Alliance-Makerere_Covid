#!/bin/bash

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
