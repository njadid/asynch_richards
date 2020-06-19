#!/bin/sh
#$ -N Richards16
#$ -j y
#$ -cwd
#$ -pe smp 6
####$ -pe 56cpn 56
####$ -l mf=16G
#$ -q IFC

/bin/echo Running on host: `hostname`.
/bin/echo In directory: `pwd`
/bin/echo Starting on: `date`

module use /Dedicated/IFC/.argon/modules
module load asynch

cd /Users/njadidoleslam/hlm_dev/asynch_richards/asynch-master
make
cd build/
make
