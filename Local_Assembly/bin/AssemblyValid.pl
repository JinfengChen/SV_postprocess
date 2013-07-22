#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"gff:s","bam:s","fa:s","step:s","project:s","help");


my $help=<<USAGE;
perl $0 --step 1
--fa : genome sequence
--gff: gff of SV to valid
  Chr1	SVpipe	Insertion	755001	755021	.	.	.	Size=-1;Method=pindel;
--bam: list of bam files for each library
  /rhome/cjinfeng/BigData/01.Rice_genomes/HEG4/00.Bam/HEG4_MSU7_BWA/FC52_7.MSU7_BWA.bam	500
--step: 1,2,3,4
--project: project dir, which will be created and store all files of this run.
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{fa}  ||="/rhome/cjinfeng/HEG4_cjinfeng/seqlib/MSU_r7.fa";
$opt{gff} ||="insertion.gff";
$opt{bam} ||="HEG4.bam.list";
$opt{project} ||= "Assembly";
`mkdir $opt{project}` unless (-e "$opt{project}");
  
my $sv=readtable($opt{gff});
my $bam=readlist($opt{bam});
my $seq=getfastaseq($opt{fa});
my $len=getfastalen($opt{fa});
my $flank=1000;
foreach my $p (sort keys %$sv){
      `mkdir $opt{project}/$p` unless (-e "$opt{project}/$p");
      my $start =$sv->{$p}->[3]-$flank >= 0 ? $sv->{$p}->[3]-$flank : 0;
      my $end   =$sv->{$p}->[4]+$flank <= $len->{$sv->{$p}->[0]} ? $sv->{$p}->[4]+$flank : $len->{$sv->{$p}->[0]};
      my $region="$sv->{$p}->[0]:$start-$end";
      my @insert=values %$bam;
      print "$p\t$region\t$start\t$end\n";
      if ($opt{step}=~/1/){ 
      ####Region sequence
      my $subseq=substr($seq->{$sv->{$p}->[0]},$start,$flank*2); 
      $subseq   =formatseq($subseq,100);
      writefile(">$region\n$subseq\n","$opt{project}/$p/Region.fa");

      ####Region reads
      foreach (@insert){
         writefile("","$opt{project}/$p/$_.raw.sam"); ### create new sam file
      }
      foreach (keys %$bam){
         `/usr/local/bin/samtools view $_ $region >> $opt{project}/$p/$bam->{$_}.raw.sam`;
      }
      foreach (@insert){
         `sort $opt{project}/$p/$_.raw.sam > $opt{project}/$p/$_.sam`;
      }
      }###step1
      
      if ($opt{step}=~/2/){
      ####Assembly
      `/rhome/cjinfeng/software/tools/Velvet/velvet/velveth $opt{project}/$p/assembly 31 -shortPaired -sam $opt{project}/$p/$insert[0].sam -shortPaired2 -sam $opt{project}/$p/$insert[1].sam`;
      `/rhome/cjinfeng/software/tools/Velvet/velvet/velvetg $opt{project}/$p/assembly -exp_cov 200 -ins_length $insert[0] -ins_length2 $insert[1] -min_contig_lgth 200 -scaffolding yes`;
      }###step2
      
      if ($opt{step}=~/3/){
      ####Alignment
      my $query="$opt{project}/$p/assembly/contigs.fa";
      my $target="$opt{project}/$p/Region.fa";
      `/opt/exonerate/2.2.0/bin/exonerate --querytype dna --targettype dna --gapextend -3 --query $query --bestn 50 --model affine:local --joinrangeext 300 --score 15 --target $target --gappedextension false --hspfilter 200 --dnahspdropoff 10 --showvulgar TRUE --showcigar TRUE --ryo "%S %pi %ql %C\n" > $opt{project}/$p/alignment.1 2> $opt{project}/$p/alignment.2`;
      }###step3

      if ($opt{step}=~/4/){
      ####Parse alignment and output refined SV


      }###step4 
}




##############################
sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/^#/);
    my @unit=split("\t",$_);
    $hash{"$unit[0]_$unit[3]_$unit[4]"}=[@unit];
}
close IN;
return \%hash;
}


sub readlist
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/^#/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}



sub writefile
{
my ($lines,$file)=@_;
open WF, ">$file" or die "$!";
     print WF "$lines";
close WF;
}


sub formatseq
{
### format a single line sequence into lines with user specific length
my ($seq,$step)=@_;
my $length=length $seq;
my $run=int ($length/$step);
my $newseq;
for(my $i=0;$i<=$run;$i++){
   my $start=$i*$step;
   my $line=substr($seq,$start,$step);
   $newseq.="$line\n";
}
return $newseq;
}



 
sub getfastaseq
{
$/=">";
my %hash;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    $hash{$head}=$seq;
}
$/="\n";
return \%hash;
}

sub getfastalen
{
$/=">";
my %hash;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    $hash{$head}= length $seq;
}
close IN;
$/="\n";
return \%hash;
}

