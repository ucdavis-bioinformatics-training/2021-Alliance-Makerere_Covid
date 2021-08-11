#!/bin/bash
#
#SBATCH --job-name=map_bwa # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=20 # Number of cores
#SBATCH --mem=16000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --time=1-00
#SBATCH --array=1
#SBATCH --partition=production # Partition to submit to
#SBATCH --account=workshop # cluster account to use for the job
#SBATCH --reservation=workshop # cluster account reservation
#SBATCH --output=slurmout/map_bowtie-%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/map_bowtie-%A_%a.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=youremail@whereever.com # Email to which notifications will be sent

start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

inpath=01-HTS_Preproc_Amplicon

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" $1`

r1=${inpath}/${sample}_R1.fastq.gz
r2=${inpath}/${sample}_R2.fastq.gz
se=${inpath}/${sample}_SE.fastq.gz

outpath='02-Bowtie2'
dname=`dirname "$sample"`
[[ -d ${outpath}/${dname} ]] || mkdir -p ${outpath}/${dname}

echo "SAMPLE: ${sample}"

THREADS=${SLURM_NTASKS}
#THREADS=20
MEM=$(expr ${SLURM_MEM_PER_NODE} / 1024)
#MEM=16
MAPTHREADS=$(expr ${THREADS} - 6)
SORTTHREADS=$(expr ${THREADS} - ${MAPTHREADS})


module load bowtie2
module load samtools

output=${outpath}/${sample}
mapfasta=./resources/NC_045512.2.fasta


## PE ##
call="bowtie2 --threads ${MAPTHREADS} \
 --very-sensitive \
 -x ${mapfasta} -1 ${r1} -2 ${r2} \
 | samtools sort -m 768M --threads ${SORTTHREADS} -o ${output}-pe.bam -"
echo $call
eval $call

## SE ##
call="bowtie2 --threads ${MAPTHREADS} \
 --very-sensitive \
 -x ${mapfasta} -U ${se} \
 | samtools sort -m 768M --threads ${SORTTHREADS} -o ${output}-se.bam -"
echo $call
eval $call

call="samtools merge -@ ${THREADS} ${output}.bam  ${output}-pe.bam ${output}-se.bam"
echo $call
eval $call

call="samtools index -@ ${THREADS} ${output}.bam"
echo $call
eval $call

call="samtools idxstats ${output}.bam > ${output}.idxstats"
echo $call
eval $call

call="samtools flagstat -@ ${THREADS} ${output}.bam > ${output}.flagstat"
echo $call
eval $call

call="samtools stats -@ ${THREADS} ${output}.bam > ${output}.stats"
echo $call
eval $call

end=`date +%s`

runtime=$((end-start))

echo $runtime
