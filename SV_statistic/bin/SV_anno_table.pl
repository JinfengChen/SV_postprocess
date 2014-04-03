#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;


GetOptions (\%opt,"list:s","project:s","help");


my $help=<<USAGE;
perl $0 --list --project
--list: list file contains GFF files of Insertion/Deletion/Inversion
Switch if to 0 first and Need to run protein_to_genome.pl first to map rice gene to SV sequence. The resulting file name as HEG4.Deletion.final.fa.genewise.gff.
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
my $output="./$opt{project}_anno";
`mkdir $output`;

my $cdslen=getfastalen("/rhome/cjinfeng/HEG4_cjinfeng/seqlib/GFF/MSU7/MSU7.all.cds");
my $cdsanno=readanno("/rhome/cjinfeng/HEG4_cjinfeng/seqlib/GFF/MSU7.gene.anno");

open IN, "$opt{list}" or die "$!";
while(<IN>){
     chomp $_;
     next if ($_=~/^$/);
     my @unit=split("\t",$_);
     my (%data,%gene,%repeat);
     $prefix=basename($unit[0],".gff");
     readsv($unit[0],\%data);
     createfasta($unit[0],"$output/$prefix.fa") unless (-e "$output/$prefix.fa");
     repeatanno("$output/$prefix.fa") unless (-e "$output/$prefix.fa.out");
     #geneanno("$output/$prefix.fa");
     if (1){
         sumrepeat("$output/$prefix.fa.out",\%data,\%repeat);
         sumgene("$output/$prefix.fa.genewise.gff",\%data,$cdslen,$cdsanno,\%gene);
         sumALL(\%repeat,\%gene,\%data,"$output/$prefix.annotation");
     }
}
close IN;

=pod
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
=cut



#########
sub repeatanno
{
my ($fasta)=@_;
`/usr/local/bin/RepeatMasker -pa 5 -species rice -nolow -no_is -norna $fasta`;
}

sub sumALL
{
my ($repeat,$gene,$all,$anno)=@_;
open ANN, ">$anno" or die "$!";
foreach my $sv (sort keys %$all){
     my $repeatanno= $repeat->{$sv} ? join("\t",@{$repeat->{$sv}}) : "NA\tNA\tNA\tNA\tNA";
     my $geneanno  = $gene->{$sv} ? join("\t",@{$gene->{$sv}}) : "NA\tNA\tNA";
     print ANN "$sv\t$all->{$sv}\t$repeatanno\t$geneanno\n"; 
}
close ANN;
}



#########
#1614   16.2  0.3  0.3  Chr10:10264822-10265112      1   291     (0) C SPMLIKE             DNA/En-Spm          (1911)   9078    8788    1
sub sumrepeat
{
my ($out,$reflen,$data)=@_;
my %hash;
open RP, "$out" or die "$!";
<RP>;
<RP>;
<RP>;
while(<RP>){
     chomp $_;
     next if ($_=~/^$/);
     my @unit=split(" ",$_);
     #print "$unit[4]\t$unit[9]\t$unit[10]\t$unit[11]\t$unit[12]\t$unit[13]\n"; 
     push @{$hash{$unit[4]}},[$unit[9],$unit[10],$unit[11],$unit[12],$unit[13]];
     
}
close RP;

print "SV\tSV_length\tTE_length\tSingle_Event\tElement:Rate\tElement:Family\n";
foreach my $sv (sort keys %hash){
     my $length=$reflen->{$sv};
     my @repeat=@{$hash{$sv}};
     my $single=0;
     if (@repeat == 1){
        my $element=shift @{$repeat[0]};
        my $family =shift @{$repeat[0]};
        my @match  =sort {$a <=> $b} @{$repeat[0]};
        my $Rleft = $1 if ($match[0]=~/\((\d+)\)/); ####left length in right end of repeat
        my $Lleft = $match[1]-1;                    ####left length in left end of repeat
        my $matchlen = $match[2]-$match[1]+1;       ####match length of repeat
        my $repeatlen= $Lleft+$matchlen+$Rleft;     ####total length of repeat
        #my $matchrate=$matchlen/$repeatlen;
        my $matchrate= sprintf("%.2f",$matchlen/$repeatlen);        ####match rate of repeat, for sv with single match > 90% means single repeat event
        my $line=join("\t",@match);
        $single = $matchrate > 0.9 ? 1 : 0;       
        print "$sv\t$length\t$matchlen\t$single\t$element:$matchrate\t$element:$family\n";
        $data->{$sv}=[$length,$matchlen,$single,"$element:$matchrate","$element:$family"];
     }elsif(@repeat == 3){
        my (@te,@temrate,@tefamily,@rate,$add);
        for(my $i=0;$i<@repeat;$i++){
           my $element=shift @{$repeat[$i]};
           push @te, $element;
           my $family =shift @{$repeat[$i]};
           my @match  =sort {$a <=> $b} @{$repeat[$i]};
           my $Rleft = $1 if ($match[0]=~/\((\d+)\)/); ####left length in right end of repeat
           my $Lleft = $match[1]-1;                    ####left length in left end of repeat
           my $matchlen = $match[2]-$match[1]+1;       ####match length of repeat
           my $repeatlen= $Lleft+$matchlen+$Rleft;     ####total length of repeat
           my $matchrate= sprintf("%.2f",$matchlen/$repeatlen);        ####match rate of repeat, for sv with single match > 90% means single repeat event
           my $line=join("\t",@match);
           push @temrate,"$element:$matchrate";
           push @tefamily,"$element:$family";
           $add+=$matchlen;
           push @rate, $matchrate;
        }
        my $mrate=join(";",@temrate);
        my $family=join(";",@tefamily);
        if ($te[0] eq $te[2] and $te[0]=~/\w+\_LTR/ and $te[1]=~/\w+\-int/ and $rate[0]+$rate[1]+$rate[2] > 2.5){
           $single=1;
        }
        print "$sv\t$length\t$add\t$single\t$mrate\t$family\n";
        $data->{$sv}=[$length,$add,$single,$mrate,$family];
     }else{
        my (@temrate,@tefamily,$add);
        for(my $i=0;$i<@repeat;$i++){
           my $element=shift @{$repeat[$i]};
           my $family =shift @{$repeat[$i]};
           my @match  =sort {$a <=> $b} @{$repeat[$i]};
           my $Rleft = $1 if ($match[0]=~/\((\d+)\)/); ####left length in right end of repeat
           my $Lleft = $match[1]-1;                    ####left length in left end of repeat
           my $matchlen = $match[2]-$match[1]+1;       ####match length of repeat
           my $repeatlen= $Lleft+$matchlen+$Rleft;     ####total length of repeat
           my $matchrate= sprintf("%.2f",$matchlen/$repeatlen);        ####match rate of repeat, for sv with single match > 90% means single repeat event
           my $line=join("\t",@match);
           push @temrate,"$element:$matchrate";
           push @tefamily,"$element:$family";
           $add+=$matchlen;
        }
        my $mrate=join(";",@temrate);
        my $family=join(";",@tefamily);
        print "$sv\t$length\t$add\t$single\t$mrate\t$family\n";
        $data->{$sv}=[$length,$add,$single,$mrate,$family];
     }
}

}



sub sumgene
{
my ($gff,$reflen,$cdslen,$cdsanno,$data)=@_;
my $refgff=parseGFFlen($gff);

print "SV\tSV_length\tTE_length\tGene:Rate\tGene:Annotation\n";
foreach my $sv (keys %$refgff){
   print ">$sv\n";
   my $svlen=$reflen->{$sv};
   my (@matchrate,@matchanno);
   foreach my $gene0 (keys %{$refgff->{$sv}}){
        my $len=$refgff->{$sv}->{$gene0};
        my $gene= $gene0=~/(.*)\-D\d+/ ? $1 : $gene0;
        my $cdslen=$cdslen->{$gene};
        my $cdsanno=$cdsanno->{$gene};
        my $matchrate=$len/$cdslen;
        $matchrate=sprintf("%.02f",$matchrate);
        push @matchrate, "$gene:$matchrate";
        push @matchanno, "$gene:$cdsanno";
        print "$gene\t$len\t$cdslen\t$cdsanno\n";
   }
   my $mrate=join(";",@matchrate);
   my $anno =join(";",@matchanno);
   print "$sv\t$svlen\t$mrate\t$anno\n";
   $data->{$sv}=[$svlen,$mrate,$anno];
}

}



#############
sub createfasta
{
my ($gff,$fasta)=@_;
open OUT,">$fasta" or die "$!";
open GFF, "$gff" or die "$!";
while(<GFF>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $id="$unit[0]:$unit[3]-$unit[4]";
    my $seq="";
    if ($unit[8]=~/Seq=(.*?);/){
       $seq=$1;
    }
    $seq=formatseq($seq,100);
    print OUT ">$id\n$seq\n";
}
close GFF;
close OUT;
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


#################
sub readsv
{
my ($gff,$data)=@_;
open GFF, "$gff" or die "$!";
while(<GFF>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $id="$unit[0]:$unit[3]-$unit[4]";
    my $len;
    if ($unit[8]=~/Seq=(.*?);/){
       $len=length $1;
    }
    $data->{$id}=$len;
}
close GFF;
}




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
    $hash{$head}=length $seq;
}
$/="\n";
return \%hash;
}



sub parseGFFlen
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;
my $id;    ##ID for element
my $record;##all line for this record
open GFF, "$gff" or die "$!";
while (<GFF>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $record="$_\n";
        $hash{$seq}{$id}=0;
    }elsif($unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$seq}{$id}+=$unit[4]-$unit[3]+1;
    }

}
close GFF;
return \%hash;
}





sub readanno
{
my ($file)=@_;
my %hash;
open AN, "$file" or die "$!";
while(<AN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=$unit[1];
}
close AN;
return \%hash;
}


 
