perl /rhome/cjinfeng/software/bin/fastaDeal.pl --sample 1 /rhome/cjinfeng/HEG4_cjinfeng/seqlib/MSU_r7.fa > Chr1.fa

perl /rhome/cjinfeng/software/bin/fastaDeal.pl -sub 754001-756021 Chr1.fa > myRegion.fa

awk '$1~/Chr1$/ || $1~/Chr10/' /rhome/cjinfeng/BigData/00.RD/Variations/SV_postprocess/SV_merge/bin/Insertion/HEG4.insertion.other.gff > insertion.gff


echo "Valid using local assembly,insertion.1"
perl AssemblyValid.pl --step 1234 > log 2> log2 &

echo "Valid using genome, insertion.0"
ln -s HEG4.insertion.other.gff insertion.gff
perl GenomeValid.pl --step 1234 > log 2> log2 &

echo "After Valid using genome, valid these need manuals using local assembly"
perl AssemblyValid.pl --gff insertion.0.Manual.gff --project Assembly.0 --step 1234 > log 2> log2 &

echo "merge Valided by genome and further valided by assembly"
cat insertion.0.LocalAssembly.gff insertion.0.Manual.1.LocalAssembly.gff | sort -k1,1 -k4,4n > insertion.ALL.LocalAssembly.gff
perl formalGFF.pl --gff insertion.ALL.LocalAssembly.gff

