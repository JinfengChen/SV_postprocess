echo "Merge Deletion"
perl SV_del_merge.pl
mv HEG4.* Deletion/
mv all.gff* Deletion/

echo "Merge Insertion"
perl SV_ins_merge.pl
mv HEG4.* Deletion/
mv all.gff* Deletion/

