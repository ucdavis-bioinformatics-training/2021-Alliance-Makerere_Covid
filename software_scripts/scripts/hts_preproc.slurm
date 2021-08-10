#!/bin/bash

#SBATCH --job-name=htstream # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --time=1-0
#SBATCH --mem=500 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # cluster partition
#SBATCH --account=workshop # cluster account to use for the job
#SBATCH --reservation=workshop # cluster account reservation
#SBATCH --array=1-24
#SBATCH --output=slurmout/htstream-%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/htstream-%A_%a.err # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@whereever.com

start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" $1`

outpath="01-HTS_Preproc_Amplicon"
dname=`dirname "$sample"`
[[ -d ${outpath}/${dname} ]] || mkdir -p ${outpath}/${dname}

echo "SAMPLE: ${sample}"

export PATH=/share/workshop/alliance_covid/software/HTStream/bin:$PATH

call="hts_Stats \
      	--notes 'initial stats' \
      	--stats-file ${outpath}/${sample}.json \
      	-1 00-RawData/${sample}_*R1*.fastq.gz \
      	-2 00-RawData/${sample}_*R2*.fastq.gz | \
      hts_Overlapper \
      	--notes 'overlap and trim adapters' \
      	--append-stats-file ${outpath}/${sample}.json | \
      hts_SeqScreener \
        --notes 'kmer match to Covid genome' \
        --append-stats-file ${outpath}/${sample}.json \
        --inverse \
      	--kmer 21 \
      	--seq resources/NC_045512.2.fasta | \
      hts_Primers \
      	--notes 'identify and trip adapters' \
      	--append-stats-file ${outpath}/${sample}.json \
      	--min_primer_matches 2 \
      	--primers_5p resources/primer5p.fa \
      	--primers_3p resources/primer3p.fa \
      	--flip | \
      hts_NTrimmer \
      	--notes 'remove any reads with  N characters' \
      	--append-stats-file ${outpath}/${sample}.json \
      	--exclude | \
      hts_LengthFilter \
      	--notes 'remove reads < 50bp' \
      	--append-stats-file ${outpath}/${sample}.json \
      	--no-orphans \
      	--min-length 50 | \
      hts_Stats \
      	--notes 'final stats' \
      	--force \
      	--append-stats-file ${outpath}/${sample}.json \
      	--fastq-output ${outpath}/${sample}"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
