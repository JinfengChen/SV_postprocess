#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;


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
        "1"   => "Whole-gene",
        "2"   => "Partial-gene",
        "3"  => "Upstream",
        "4"  => "Utr5prime",
        "5" => "Exon",
        "6" => "Intron",
        "7"=> "Utr3prime",
        "8"=> "Downstream"
);

$opt{project} ||= "HEG4";
my $output="./$opt{project}_anno";
`mkdir $output` unless (-e "$output");
my $refanno=readanno("/rhome/cjinfeng/HEG4_cjinfeng/seqlib/GFF/MSU7.gene.anno");
my $refinf =readinf("/rhome/cjinfeng/HEG4_cjinfeng/seqlib/GFF/RiceGene/MSU7.gene.inf");

my %data;
my @type;
my $BIN="/rhome/cjinfeng/software/tools/SVcaller/snpEff/";
open IN, "$opt{list}" or die "$!";
while(<IN>){
     chomp $_;
     next if ($_=~/^$/);
     my @unit=split("\t",$_);
     my $prefix=basename($unit[0],".gff");
     push @type, $prefix;
     GFF2BED($unit[0],"$output/$prefix.bed");
     `/opt/java/jdk1.6.0_38/bin/java -Xmx1g -jar $BIN/snpEff.jar eff -c $BIN/snpEff.config -v -i bed -o bed rice7 $output/$prefix.bed > $output/$prefix.eff.bed` unless (-e "$output/$prefix.eff.bed");
     readeffbed("$output/$prefix.eff.bed", "$output/$prefix.eff.sum",$prefix,\%data,$refanno,$refinf);
}
close IN;

open OUT, ">$opt{project}.SV.insertion_site.distr" or die "$!";
print OUT "Effect\t$type[0]\t$type[1]\n";
foreach my $rank (sort keys %hash){
     my $t=$hash{$rank};
     print OUT "$t";
     for(my $i=0;$i<@type;$i++){
        my $num=$data{$t}{$type[$i]} ? $data{$t}{$type[$i]} : 0;
        print OUT "\t$num";
     }
     print OUT "\n";
}
close OUT;

draw();

#`rm $output/$opt{project}.SV.bed`;
#`cat $output/*.bed > $output/$opt{project}.SV.bed`;
#my $BIN="/rhome/cjinfeng/software/tools/SVcaller/snpEff/";
#`/opt/java/jdk1.6.0_38/bin/java -Xmx1g -jar $BIN/snpEff.jar eff -c $BIN/snpEff.config -v -i bed -o bed rice7 $output/$opt{project}.SV.bed > $output/$opt{project}.SV.eff.bed` unless (-e "$output/$opt{project}.SV.eff.bed");



sub readanno
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $unit[0] = $1 if ($unit[0]=~/(.*?)\.\d+/);
    $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}

sub readinf
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[2]}=$_;
}
close IN;
return \%hash;
}




######################################
#1       100145  103713  Deletion_1;Upstream|LOC_Os01g01210.1||LOC_Os01g01210|mRNA;Exon|exon_1_2|LOC_Os01g01210.1||LOC_Os01g01210|mRNA;
sub readeffbed
{
my ($bed,$sum,$prefix,$data,$refanno,$refinf)=@_;
my %hash;
open BED, "$bed" or die "$!";
while(<BED>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/^#/);
    my @unit=split("\t",$_);
    my @anno=split(";",$unit[3]);
    for(my $i=1;$i<@anno;$i++){
       my ($eff,$mrna);
       if($anno[$i]=~/(.*?)\|.*?\|(.*?)\|\|(.*?)\|.*?$/){ ###Exon|exon_1_2|LOC_Os01g01210.1||LOC_Os01g01210|mRNA
          $eff=$1;
          $mrna=$3;
          $hash{$anno[0]}{$mrna}{$eff}=1;
       }elsif ($anno[$i]=~/(.*?)\|(.*?)\|\|(.*?)\|.*?$/){ ###Upstream|LOC_Os01g01210.1||LOC_Os01g01210|mRNA
          $eff=$1;
          $mrna=$3;
          $hash{$anno[0]}{$mrna}{$eff}=1;
       }
    }
}
close BED;
open SUM, ">$sum" or die "$!";
foreach my $sv (sort keys %hash){ 
foreach my $gene (sort keys %{$hash{$sv}}){
    next unless (exists $refanno->{$gene});
    my $inf = $refinf->{$gene} ? $refinf->{$gene} : "NA";
    my @eff=keys %{$hash{$sv}{$gene}};
    my $effect=join(";",@eff);
    if (@eff > 1){
       if ($effect=~/Upstream/ and $effect=~/Exon/ and  $effect=~/Downstream/){
          print SUM "$sv\t$gene\tWhole-gene\t$effect\t$refanno->{$gene}\t$inf\n";
          $data{"Whole-gene"}{$prefix}++;
       }else{
          print SUM "$sv\t$gene\tPartial-gene\t$effect\t$refanno->{$gene}\t$inf\n";
          $data{"Partial-gene"}{$prefix}++;
       }
    }else{
       print SUM "$sv\t$gene\t$eff[0]\t$eff[0]\t$refanno->{$gene}\t$inf\n";
       $data{$eff[0]}{$prefix}++;
    }
}
}
close SUM;
}



#######################################
sub writefile
{
my ($lines,$file)=@_;
open WF, ">$file" or die "$!";
     print WF "$lines";
close WF;
}

#######################################
#Chr1	SVpipe	Deletion	100145	103713	.	.	.	Size
sub GFF2BED
{
my ($gff,$bed)=@_;
my $count=0;
open OUT, ">$bed" or die "$!";
open GFF, "$gff" or die "$!";
while(<GFF>){
    chomp $_;
    next if ($_=~/^$/);
    $count++;
    my @unit=split("\t",$_);
    my $name=$unit[2]."_$count";
    print OUT "$unit[0]\t$unit[3]\t$unit[4]\t$name\t.\t+\n";
}
close GFF;
close OUT;
}


###Chr10:10575767-10576566 800     800     797     1       CACTA-G1:1.00   CACTA-G1:DNA/En-Spm     NA      NA      NA

sub singleTE
{
my ($data,$prefix,$anno)=@_;
open AN, "$anno" or die "$!";
while(<AN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    if ($unit[4] == 1){
       my @family=split(";",$unit[6]);
       #print "$_\n$family[0]\n";
       my $tef=$1 if ($family[0]=~/.*?\:(.*?)$/);
       $data->{$tef}->{$prefix}++;
    }
}
close AN;
}

sub draw
{
my $cmd =<<R;
pdf("$opt{project}.SV.insertion_site.pdf")
par(mar=c(10,4,4,2))
x <- read.table("$opt{project}.SV.insertion_site.distr",skip=1)
data <- rbind(x[,2]/sum(x[,2]),x[,3]/sum(x[,3]))
xx <- barplot(data,beside=TRUE,ylab="Proportion",border=FALSE,ylim=c(0,0.5),col=c("Orange","blue"))
axis(1,c(0.5,max(xx)+0.5),line=0,labels=c("",""))

for (i in 1:length(xx[1,])) { # adj can not take vector, so we use loops to add text
  text(xx[1,i],-0.03,labels=x[i,1],cex=1,srt=55,adj=c(1,1),xpd=TRUE)
}

text(xx[1,],x[,2]/sum(x[,2])+0.03,offset=2,labels=x[,2],srt=55,xpd=TRUE)
text(xx[1,]+1.2,x[,3]/sum(x[,3])+0.03,offset=2,labels=x[,3],srt=55,xpd=TRUE)
legend("topleft",c("Deletion","Insertion"),bty="n",border="NA",lty=c(0,0),cex=1,fill=c("orange","blue"))
mtext("Effect on Gene",side=1, at=12,line=8)
dev.off()

R

open OUT, ">$opt{project}.SV.insertion_site.R" or die "$!";
     print OUT "$cmd\n";
close OUT;
`cat $opt{project}.SV.insertion_site.R | R --vanilla --slave`;
}


