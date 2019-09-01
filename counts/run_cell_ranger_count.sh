#!/bin/bash

set -v
set -e 

echo $PWD

module load bcl2fastq/2.18
module load cellranger

fastqfolder=../fastq/
transcriptome=/nfs/rprdata/refGenome10x/refdata-cellranger-hg19-1.2.0/

samplefile=$fastqfolder/outs/input_samplesheet.csv

find ${fastqfolder} -name '*fastq.gz' | sed 's/.*\///;s/_S1.*//' | sort | uniq > libList.txt


## Submit cell-ranger jobs to the grid. 
cat libList.txt | \
while read sample; 
do 
    echo "#################"
    echo $sample 
    echo $fastqlist
    echo "cd $PWD; module load bcl2fastq cellranger; 
time cellranger count \
      --id=$sample \
      --fastqs=$fastqfolder \
      --sample=$sample \
      --transcriptome=$transcriptome \
      --localcores=24 --localmem=110" | qsub -q erprq -l nodes=1:ppn=28 -l mem=120g -N $sample
    sleep 0.5;
done


