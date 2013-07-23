#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=5gb
#PBS -l walltime=10:00:00

cd $PBS_O_WORKDIR

#perl AssemblyValid.pl --step 1234
perl GenomeValid.pl --step 1234

echo "Done"

