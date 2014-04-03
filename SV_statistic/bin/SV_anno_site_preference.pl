#!/usr/bin/perl
use Getopt::Long;
use strict;
use warnings;
use FindBin qw($Bin);

my %opt;

GetOptions (\%opt,"gene:s","repeat:s","utr:s","trf:s","unalign:s","project:s","help");


my $help=<<USAGE;
perl $0 --gene gene.gff --repeat repeat.gff --utr utr.gff --trf trf.gff --unalign unalign.table --project BAC
For insertion this probably a preference site. But for deletion, this could be a influence on genome.
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $te=parseGFFrepeat($opt{repeat});
& gff2table($opt{gene},"CDS","$opt{project}.CDS");
& gff2table($opt{gene},"intron","$opt{project}.intron");
& gff2table($opt{utr},"UTR","$opt{project}.utr");
& gff2table($opt{trf},"TandemRepeat","$opt{project}.trf");
& gff2table($opt{repeat},"Trans","$opt{project}.repeat");

`perl $Bin/scripts/findOverlap.pl $opt{unalign} $opt{project}.CDS > $opt{project}.CDS.overlap`;
`perl $Bin/scripts/findOverlap.pl $opt{unalign} $opt{project}.intron > $opt{project}.intron.overlap`;
`perl $Bin/scripts/findOverlap.pl $opt{unalign} $opt{project}.utr > $opt{project}.utr.overlap`;
`perl $Bin/scripts/findOverlap.pl $opt{unalign} $opt{project}.trf > $opt{project}.trf.overlap`;
`perl $Bin/scripts/findOverlap.pl $opt{unalign} $opt{project}.repeat > $opt{project}.repeat.overlap`;


my ($cds)=sumCDS("$opt{project}.CDS.overlap");
my ($intron)=sumCDS("$opt{project}.intron.overlap");
my ($utr)=sumCDS("$opt{project}.utr.overlap");
my ($trf)=sumCDS("$opt{project}.trf.overlap");
my $repeat=sumCDS("$opt{project}.repeat.overlap",$te);

print "$cds\t$intron\t$utr\t$trf\t$repeat\n";

=cut
print "Total Length:$total\n";
my $cdsrate=$cds/$total;
print "CDS: $cds\t$cdsrate\n";
my $intronrate=$intron/$total;
print "Intron: $intron\t$intronrate\n";
my $utrrate=$utr/$total;
print "UTR: $utr\t$utrrate\n";
my $trfrate=$trf/$total;
print "TRF: $trf\t$trfrate\n";
my $junk;
foreach my $t (sort keys %$repeat){
   $junk+=$repeat->{$t};
   my $rate=$repeat->{$t}/$total;
   print "$t\t$repeat->{$t}\t$rate\n";
}
my $unknown=$total-$junk-$cds-$intron-$utr-$trf;
my $unknownrate=$unknown/$total;
print "Unknown: $unknown\t$unknownrate\n";
=cut

###################
sub parseGFFrepeat
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[8]=~/ID=(.*?);.*Class=(.*?);/){
       my $id=$1;
       my $class=$2;
       if ($class=~/LTR/){
          $class="LTR";
       }elsif($class=~/DNA/){
          $class="DNA";
       }else{
          $class="otherTE";
       }
       $hash{$id}=$class
    }
}
close IN;
return \%hash;
}


#HEG4.SV.Site_Preference.Deletion.CDS.overlap
#Deletion_18     1317    Chr1    2       LOC_Os01g05380.1,27,27  LOC_Os01g05380.1,210,210
sub sumCDS
{
my ($file)=@_;
my $count =0 ;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split(" ",$_);
    if ($unit[3] > 0){
       $count++;
    }
}
close IN;
return ($count);
}

sub sumintron
{
my ($file)=@_;
my $incdslen;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split(" ",$_);
    next if ($unit[3] == 0);
    for(my $i=4;$i<@unit;$i++){
       my @array=split(",",$unit[$i]);
       $incdslen+=$array[2];
    }
}
close IN;
return ($incdslen);
}



sub sumrepeat
{
my ($file,$te)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split(" ",$_);
    next if ($unit[3] == 0);
    my %temp;
    for(my $i=4;$i<@unit;$i++){
       my @array=split(",",$unit[$i]);
       my $id=$array[0];
       my $type=$te->{$id};
       my $len=$array[2];
       $hash{$type}+=$len;
    }
}
close IN;
return \%hash;
}

sub gff2table
{
my ($gff,$feature,$table)=@_;
open IN, "$gff" or die "$!";
open OUT, ">$table.unsort" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "" or $_ =~ /^#/);
    my @unit=split("\t",$_);
    if ($unit[2] =~ $feature and ($unit[8] =~/ID=(.*?);/ or $unit[8] =~/Parent=(.*)/)){
        #$unit[0]=~tr/[A-Z]/[a-z]/;
        print OUT "$unit[0]\t$1\t$unit[3]\t$unit[4]\n";
    }
}
close IN;
close OUT;
system ("msort -k 1,n3 $table.unsort > $table");
system ("rm $table.unsort");
}

