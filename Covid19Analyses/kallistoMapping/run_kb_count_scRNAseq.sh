#!/bin/bash

set -v
set -e 

echo $PWD


## Can use this one as an argument. 
##fastqfolder=../../fastq/
fastqfolder=${PWD}/../fastq/

kref=${PWD}/ref/

find ${fastqfolder} -name 'HPL*fastq.gz' | sed 's/.*\///;s/_S.*//' | sort | uniq > libList.txt

mkdir -p bus

cat libList.txt | \
while read sample; 
do 
##    fastqs=`find ${fastqfolder} -name "${sample}*R[12]*fastq.gz"`
##    fastqlist=`echo ${fastqs} | tr ' ' ,`
##time kallisto bus -i ${kallistoref} -o ./bus/${sample} -x 10xv2 -t 20 ${fastqs}
    echo "#################"
    fastqs=`find ${fastqfolder} -name "${sample}*R[12]*fastq.gz" | sort | uniq`
    fastqlist=`echo ${fastqs} | tr ' ' ' '`
    echo "#################"
    echo $sample 
    echo $fastqlist
##    echo $fastqs
    echo "cd $PWD; 
module unload python;
module load anaconda3.python;
source activate kallisto2;
mkdir -p \${TMPDIR}/${sample};
cd \${TMPDIR}/${sample};
time kb count --h5ad -i ${kref}/index.idx -o ${sample} -g ${kref}/t2g.txt -x 10xv3 -c1 ${kref}/spliced_t2c.txt -c2 ${kref}/unspliced_t2c.txt --lamanno --filter bustools --verbose -t 6 ${fastqlist};
mv ${sample} ${PWD}/bus/;" | qsub -q erprq -l nodes=1:ppn=8 -l mem=90g -N $sample
    sleep 0.5;
done


# kb count --h5ad -i index.idx -g t2g.txt -x 10xv2 -o SRR6470906 \
# -c1 spliced_t2c.txt -c2 unspliced_t2c.txt --lamanno --filter bustools -t 2 \
# SRR6470906_S1_L001_R1_001.fastq.gz \
# SRR6470906_S1_L001_R2_001.fastq.gz \
# SRR6470906_S1_L002_R1_001.fastq.gz \
# SRR6470906_S1_L002_R2_001.fastq.gz


## --jobmode=erprq
## qsub -I -q erprq -l nodes=1:ppn=28 -l mem=120g
