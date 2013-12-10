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
#my $refanno=readanno("/rhome/cjinfeng/HEG4_cjinfeng/seqlib/GFF/MSU7.gene.anno");
#my $refinf =readinf("/rhome/cjinfeng/HEG4_cjinfeng/seqlib/GFF/RiceGene/MSU7.gene.inf");
my $reffpkm =readfpkm("/rhome/cjinfeng/BigData/02.Transcription/Transcriptome/bin/cufflink.table");

my %data;
my @type;
open IN, "$opt{list}" or die "$!";
while(<IN>){
     chomp $_;
     next if ($_=~/^$/);
     my @unit=split("\t",$_);
     my $prefix=basename($unit[0],".gff");
     push @type, $prefix; 
     readeffsum("$output/$prefix.eff.sum",$prefix,\%data,$reffpkm);
}
close IN;

open OUT, ">$opt{project}.SV.FPKM.distr" or die "$!";
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

#draw();

#`rm $output/$opt{project}.SV.bed`;
#`cat $output/*.bed > $output/$opt{project}.SV.bed`;
#my $BIN="/rhome/cjinfeng/software/tools/SVcaller/snpEff/";
#`/opt/java/jdk1.6.0_38/bin/java -Xmx1g -jar $BIN/snpEff.jar eff -c $BIN/snpEff.config -v -i bed -o bed rice7 $output/$opt{project}.SV.bed > $output/$opt{project}.SV.eff.bed` unless (-e "$output/$opt{project}.SV.eff.bed");



sub readfpkm
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $unit[0] = $1 if ($unit[0]=~/(.*?)\.\d+/);
    my @fpkm=($unit[2],$unit[3]);
    my $sd=sd(\@fpkm); 
    #print "$unit[0]\t$sd\n";
    $hash{$unit[0]}=$sd;
}
close IN;
return \%hash;
}

sub sd
{
my ($num)=@_;
my $loop=0;
my $total=0;
my $add_square;
foreach  (@$num) {
        next if ($_ eq "NA");
        #my $temp=log($_);
        my $temp=$_;
        $total+=$temp;
        $add_square+=$temp*$temp;
        $loop++;
}
my $number=$loop;
return (0) if ($number < 2);
my $mean=$total/$number;
my $SD=sqrt( ($add_square-$total*$total/$number)/ ($number-1) );
my $se=1.96*$SD/sqrt($number);
return ($SD);
}

sub mean
{
my ($num)=@_;
my $loop=0;
my $total;
foreach  (@$num) {
        next if ($_ eq "NA");
        $total+=$_;
        $loop++;
}
return 0 if ($loop == 0);
my $mean=$total/$loop;
return $mean;
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
#Deletion_10     LOC_Os01g04660  Downstream      Downstream      lipid phosphatase protein,putative,expressed    NA
sub readeffsum
{
my ($sum,$prefix,$data,$reffpkm)=@_;
my %hash;
open BED, "$sum" or die "$!";
while(<BED>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/^#/);
    my @unit=split("\t",$_);
    push @{$hash{$prefix}{$unit[2]}},$reffpkm->{$unit[1]};    
}
close BED;

foreach my $pre (keys %hash){
    foreach my $type (keys %{$hash{$pre}}){
        my @num=@{$hash{$pre}{$type}};
        my $mean=mean(\@num);
        print "$pre\t$type\t$mean\n";
        $data->{$type}->{$pre}=$mean;
    }
}

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


