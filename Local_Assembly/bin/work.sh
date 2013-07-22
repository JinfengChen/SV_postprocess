perl /rhome/cjinfeng/software/bin/fastaDeal.pl --sample 1 /rhome/cjinfeng/HEG4_cjinfeng/seqlib/MSU_r7.fa > Chr1.fa

perl /rhome/cjinfeng/software/bin/fastaDeal.pl -sub 754001-756021 Chr1.fa > myRegion.fa

awk '$1~/Chr1$/ || $1~/Chr10/' /rhome/cjinfeng/BigData/00.RD/Variations/SV_postprocess/SV_merge/bin/Insertion/HEG4.insertion.other.gff > insertion.gff

perl AssemblyValid.pl --step 1234 > log 2> log2 &


