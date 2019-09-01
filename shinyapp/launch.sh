#!/bin/bash
## Script is submitted to this Queue:
#PBS -q wsuq

#PBS -l select=1:ncpus=4:mem=10GB

set -e 
set -x

module load R

cd /nfs/rprdata/scilab/novogene/Analyses/ED6/PrimaryAnalysis/shinyapp_Roger

hostIP=`hostname -i`

piquelabPort=9091

## Preparing tunnel 

ssh piquelab ssh -g -4 -nTL ${piquelabPort}:${hostIP}:4236 -N warrior &

echo "http://piquelab.grid.wayne.edu:"${piquelabPort}

echo "shiny::runApp('./',port=4236,host='"${hostIP}"')" | R --vanilla 


