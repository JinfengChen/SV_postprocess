#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"help");


my $help=<<USAGE;
perl $0
Merge insertion predicted by pindel, bd, meerkat, soapsv. Skip these >= 100 kb because of high false positive by manual inspection.
Cluster insertions by BEDtools (cluster) and selected representative coordinate ranked by soapsv, pindel, meerkat and bd.
Pindel and soapsv are basepair resolution and meerkat using insertion size and split reads also should be basepair resolution. bd use only insertion size.
USAGE


if ($opt{help}){
    print "$help\n";
    exit();
}

my @gff=glob("../input/Insertion/HEG4.*.gff");
my @del;
for(my $i=0;$i<@gff;$i++){
   unless ($gff[$i]=~/mping/i){
      my $prefix=basename($gff[$i],".gff");
      `grep "Insertion" $gff[$i] | awk '{len=\$5-\$4;if(len < 100000){print}}' > $prefix.insertion.gff` unless (-e "$prefix.insertion.gff");
      push @del, "$prefix.insertion.gff";
   }else{
      `cat $gff[$i] | sed 's/HEG4_2\.3\.mPing\.20X/RelocaTE/' | sed 's/transposable_element_insertion_site/Insertion/' | sed 's/ID=.*/Size=430;/' > HEG4.mping.insertion.gff`;
      push @del, "HEG4.mping.insertion.gff";
   }
}
#my $deletion=join(" ", @del);
#`bash /rhome/cjinfeng/HEG4_cjinfeng/Variations/VCFutility/bin/mergeGFF.sh HEG4.Deletion.gff $deletion`;

#Draw overlap among different method, all.gff will be used to cluster deletions and merge them.
my $deletion=join(",",@del);
`perl /rhome/cjinfeng/HEG4_cjinfeng/Variations/VCFutility/bin/GFF2venn.pl --gff $deletion` unless (-e "all.gff");

#Cluster deletions
`/opt/bedtools/2.17.0-25-g7b42b3b/bin//bedtools cluster -i all.gff > all.gff.cluster` unless (-e "all.gff.cluster");
`cat *.insertion.gff > HEG4.all.gff` unless (-e "HEG4.all.gff");
my $cluster=readcluster("all.gff.cluster");
my $gff=readgff("HEG4.all.gff");

#merge deletions
my $title="SVpipe";
open OUT1, ">HEG4.insertion.known.gff" or die "$!";  ###gff with insertion sequence known by soapsv or TE tools like RelocaTE
open OUT2, ">HEG4.insertion.other.gff" or die "$!";  ###gff with insertion sequence unkonwn, need to do local assembly
foreach my $c (sort {$a <=> $b} keys %$cluster){
   my @d=@{$cluster->{$c}};
   my %hash;
   for(my $i=0;$i<@d;$i++){ ### for each deletion in each cluster
      if ($d[$i]->[3]=~/pindel/i){
         $hash{"pindel"}=[$d[$i]->[3],"$d[$i]->[0]\_$d[$i]->[1]\_$d[$i]->[2]"];
      }elsif($d[$i]->[3]=~/soapsv/i){
         $hash{"soapsv"}=[$d[$i]->[3],"$d[$i]->[0]\_$d[$i]->[1]\_$d[$i]->[2]"];
      }elsif($d[$i]->[3]=~/meerkat/i){
         $hash{"meerkat"}=[$d[$i]->[3],"$d[$i]->[0]\_$d[$i]->[1]\_$d[$i]->[2]"];
      }elsif($d[$i]->[3]=~/breakdancer/i){
         $hash{"breakdancer"}=[$d[$i]->[3],"$d[$i]->[0]\_$d[$i]->[1]\_$d[$i]->[2]"];
      }elsif($d[$i]->[3]=~/relocate/i){
         $hash{"relocate"}=[$d[$i]->[3],"$d[$i]->[0]\_$d[$i]->[1]\_$d[$i]->[2]"];
      }else{
         $hash{"new"}=[$d[$i]->[3],"$d[$i]->[0]\_$d[$i]->[1]\_$d[$i]->[2]"];
      }
   }
   my $method;
   if (exists $hash{"relocate"}){
       $method="relocate";
   }elsif(exists $hash{"soapsv"}){
      #print OUT "$c\t$hash{pindel}->[0]\t$hash{pindel}->[1]\n";
      #print OUT "$c\t$gff->{$hash{pindel}->[0]}{$hash{pindel}->[1]}\n";
       $method="soapsv";
   }elsif(exists $hash{"pindel"}){
      #print OUT "$c\t$hash{soapsv}->[0]\t$hash{soapsv}->[1]\n";
      #print OUT "$c\t$gff->{$hash{soapsv}->[0]}{$hash{soapsv}->[1]}\n";
       $method="pindel";
   }elsif(exists $hash{"meerkat"}){
      #print OUT "$c\t$hash{meerkat}->[0]\t$hash{meerkat}->[1]\n";
      #print OUT "$c\t$gff->{$hash{meerkat}->[0]}{$hash{meerkat}->[1]}\n";
       $method="meerkat";
   }elsif(exists $hash{"breakdancer"}){
      #print OUT "$c\t$hash{breakdancer}->[0]\t$hash{breakdancer}->[1]\n";
      #print OUT "$c\t$gff->{$hash{breakdancer}->[0]}{$hash{breakdancer}->[1]}\n";
       $method="breakdancer";
   }else{
       $method="new";
   }
   my @temp=split("\t",$gff->{$hash{$method}->[0]}{$hash{$method}->[1]});
   $temp[1]=$title;
   $temp[8].="Method=".$method.";";
   my $g=join("\t",@temp);
   if ($method=~/soapsv/ or $method=~/relocate/){
      print OUT1 "$g\n";
   }else{
      print OUT2 "$g\n";
   }
}
close OUT1;
close OUT2;

################
#Chr1    Meerkat2        Deletion        100147  103717  .       .       .       Size=3569;Mech=TEI_LTR/Gypsy_SZ-40_LTR;Type=del;
sub readgff
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/^#/);
    my @unit=split("\t",$_);
    $hash{$unit[1]}{"$unit[0]\_$unit[3]\_$unit[4]"}=$_;
}
close IN;
return \%hash;
}

#Chr1    100145  103713  Pindel1 1
sub readcluster
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    push @{$hash{$unit[4]}}, [@unit];
}
close IN;
return \%hash;
}

