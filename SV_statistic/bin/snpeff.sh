#!/bin/sh

###Annotation
REFERENCE="rice7"
DB=/rhome/cjinfeng/HEG4_cjinfeng/Variations/SV/SNPeff/input/
BIN=/rhome/cjinfeng/software/tools/SVcaller/snpEff
INPUT_BED=./HEG4_anno/HEG4.SV.bed
BED_BASE=`basename $INPUT_BED .bed`
SNPEFF="/opt/java/jdk1.6.0_38/bin/java -Xmx1g -jar $BIN/snpEff.jar"

$SNPEFF eff -c $BIN/snpEff.config -v -i bed -o bed $REFERENCE $INPUT_BED > ./HEG4_anno/$BED_BASE.eff.bed


