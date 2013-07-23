#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"gff:s","help");


my $help=<<USAGE;
perl $0 --gff 
Format information feild of gff to Size=#;Method=***;INDEL=Insertion;Seq=***;
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

readtable("$opt{gff}");


##Chr1    SVpipe  Insertion       2566394 2566396 .       .       .       Size=-1;Method=pindel;INDEL=Insertion;Seq=
sub readtable
{
my ($file)=@_;
my $prefix=basename($file,".gff");
open OUT, ">$prefix.final.gff" or die "$!";
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $method="NA";
    if ($unit[8]=~/Method=(.*?);/){
       $method=$1;
    }
    my $indel="NA";
    if ($unit[8]=~/INDEL=(.*?);/){
       $indel=$1;
    }
    if ($unit[8]=~/Seq=(.*?);/){
       my $seq=$1;
       my $size=length $1;
       $unit[8]="Size=$size;Method=$method;INDEL=$indel;Seq=$seq;";
    }
    my $line=join("\t",@unit);
    print OUT "$line\n";
}
close IN;
close OUT;
}
 
