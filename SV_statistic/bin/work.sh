echo "create SV file list"
ls `pwd`/../input/HEG4* > HEG4.SV.list
echo "length distribution"
perl SV_len_distr.pl --list HEG4.SV.list
cat length.R | R --slave

