#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"list:s","project:s","help");


my $help=<<USAGE;
perl $0 --list --project
--list: list file contains GFF files of Insertion/Deletion/Inversion
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my %hash=(
        "1"   => "1",
        "2"   => "2-10",
        "10"  => "10-50",
        "50"  => "50-100",
        "100" => "100-500",
        "500" => "500-1000",
        "1000"=> "1000-5000",
        "5000"=> ">5000"
);

$opt{project} ||= "HEG4";

my %data;
open IN, "$opt{list}" or die "$!";
while(<IN>){
     chomp $_;
     next if ($_=~/^$/);
     my @unit=split("\t",$_);
     distr($unit[0],\%data);
}
close IN;

#print "Rank\tRange\tALL\tRefInsertion\tQryInsertion\n";
#foreach my $k (sort {$a <=> $b} keys %$all){
#   print "$k\t$hash{$k}\t$all->{$k}\t$rice->{$k}\t$oga->{$k}\n";
#}

open OUT, ">$opt{project}.SV.length.distr" or die "$!";
print OUT "Rank\tRange";
foreach my $type (sort keys %data){
      print OUT "\t$type";
}
print OUT "\n";

foreach my $rank (sort {$a <=> $b} keys %hash){
     print OUT "$rank\t$hash{$rank}";  
     foreach my $type (sort keys %data){
           my $num=$data{$type}{$rank} > 0 ? $data{$type}{$rank} : 0;
           print OUT "\t$num";
     }
     print OUT "\n";
}
close OUT;




##############
sub distr
{
my ($file,$hash)=@_;
open DI, "$file" or die "$!";
while(<DI>){
     chomp $_;
     next if ($_=~/^$/);
     my @unit=split("\t",$_);
     my $len = $1 if ($unit[8]=~/Size=(\d+?);/);
     if ($len == 1){
        $hash->{$unit[2]}->{"1"}++;
     }elsif($len < 10){
        $hash->{$unit[2]}->{"2"}++;
     }elsif($len < 50){
        $hash->{$unit[2]}->{"10"}++;
     }elsif($len < 100){
        $hash->{$unit[2]}->{"50"}++;
     }elsif($len < 500){
        $hash->{$unit[2]}->{"100"}++;
     }elsif($len < 1000){
        $hash->{$unit[2]}->{"500"}++;
     }elsif($len < 5000){
        $hash->{$unit[2]}->{"1000"}++;
     }else{
        $hash->{$unit[2]}->{"5000"}++;
     }
}
close DI;
} 
