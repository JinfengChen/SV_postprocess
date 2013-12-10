echo "create SV file list"
ls `pwd`/../input/HEG4* > HEG4.SV.list

echo "length distribution"
perl SV_len_distr.pl --list HEG4.SV.list
cat length.R | R --slave

echo "annotation"
perl SV_anno_table.pl --list HEG4.SV.list > log 2> log2 &

echo "after SV_anno_table.pl, run single TE insertion diff between NB and HEG4"
perl SV_anno_TEinsertion.pl --list HEG4.SV.list

echo "annotation insertion site in reference"
perl SV_anno_insertion_site.pl --list HEG4.SV.list

