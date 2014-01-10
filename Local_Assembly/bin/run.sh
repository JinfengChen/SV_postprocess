#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=12gb
#PBS -l walltime=10:00:00

cd $PBS_O_WORKDIR

echo "Deletion"
#perl GenomeValid_Deletion.pl --gff insertion.gff --step 1234
#perl AssemblyValid_Deletion.pl --gff insertion.0.Manual.gff --project Assembly.0 --step 1234

echo "Insertion"
perl GenomeValid.pl --step 1234
perl AssemblyValid.pl --gff insertion.0.Manual.gff --project Assembly.0 --step 1234

echo "Done"

