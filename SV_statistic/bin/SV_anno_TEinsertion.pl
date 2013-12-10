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

my %data;
my @type;
open IN, "$opt{list}" or die "$!";
while(<IN>){
     chomp $_;
     next if ($_=~/^$/);
     my @unit=split("\t",$_);
     $prefix=basename($unit[0],".gff");
     singleTE(\%data,$prefix,"$output/$prefix.annotation");
     push @type,$prefix;
}
close IN;

open OUT, ">$opt{project}.SV.singleTE.distr" or die "$!";
print OUT "TE_family\t$type[0]\t$type[1]\n";
foreach my $te (sort keys %data){
     print OUT "$te";
     for(my $i=0;$i<@type;$i++){ 
        my $temp = $data{$te}{$type[$i]} ? $data{$te}{$type[$i]} : "0";
        print OUT "\t$temp";
     }
     print OUT "\n";
}
close OUT;

draw();



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
pdf("$opt{project}.SV.singleTE.pdf")
par(mar=c(10,4,4,2))
x <- read.table("HEG4.SV.singleTE.distr",skip=1)
data <- rbind(x[,2]/sum(x[,2]),x[,3]/sum(x[,3]))
xx <- barplot(data,beside=TRUE,ylab="Proportion",border=FALSE,ylim=c(0,1),col=c("Orange","blue"))
axis(1,c(0.5,max(xx)+0.5),line=0,labels=c("",""))

for (i in 1:length(xx[1,])) { # adj can not take vector, so we use loops to add text
  text(xx[1,i],-0.06,labels=x[i,1],cex=1,srt=55,adj=c(1,1),xpd=TRUE)
}

text(xx[1,],x[,2]/sum(x[,2])+0.03,offset=2,labels=x[,2],srt=55,xpd=TRUE)
text(xx[1,]+1.2,x[,3]/sum(x[,3])+0.03,offset=2,labels=x[,3],srt=55,xpd=TRUE)
legend("topright",c("Deletion","Insertion"),bty="n",border="NA",lty=c(0,0),cex=1,fill=c("orange","blue"))
mtext("Repeat Family",side=1, at=24,line=8)
dev.off()

R

open OUT, ">$opt{project}.SV.singleTE.R" or die "$!";
     print OUT "$cmd\n";
close OUT;
`cat $opt{project}.SV.singleTE.R | R --vanilla --slave`;
}


