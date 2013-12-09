

my ($indel,$sequence)=parse_age("alignment.2","Region.fa","contigs.fa");

print "$indel,$sequence";

###
#Alignment:
# first  seq =>  [  2,  9] EXCISED REGION [ 10,173]
# second seq =>  [174,168] EXCISED REGION [165,  2]
sub parse_age
{
my ($align,$target,$query)=@_;

my $indel=3;
my $sequence="NA";

my $refseq=getfastaseq($query);
$/="\n\n\n";
open IN, "$align" or die "$!";
while(<IN>){
   my $block=$_;
   while($block=~/(MATCH.*Alignment time is .*? s)/sg){  ### match every alignment block for each contig
      #print "MATCH\n$1\nEND\n";
      my $match=$1;
      my $contig;
      if ($match=~/Second seq .* nucs \'(.*)\'\n/){  ### match contig name
         $contig=$1;
         print "$contig\n";
      }
      ### no excise
      if ($match=~/Alignment:\n first\s+seq =>\s+\[\s*(\d+)\,\s*(\d+)\]\n\s+second seq \=\> \s+\[\s*(\d+)\,\s*(\d+)\]\n/){
         print "$1\t$2\t$3\t$4\n";
         if ($1 < 1000-100 and $2 > 1000+100){ ### perfect alignment cover 200 bp of breakpoint
            $indel=0 unless ($indel == 1); ### indicate no SV exists
         }
      ### excise
      }elsif($match=~/Alignment:\n first\s+seq =>\s+\[\s*(\d+)\,\s*(\d+)\] EXCISED REGION \[\s*(\d+)\,\s*(\d+)\]\n\s+second seq \=\> \s+\[\s*(\d+)\,\s*(\d+)\] EXCISED REGION \[\s*(\d+)\,\s*(\d+)\]\n/){
         print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\n";
         if ($1 < 1000-100 and $4 > 1000+100){ ### excised alignment cover 200 bp of breakpoint
            if ($2 <= 1000+100 and $2 >= 1000-100 and $3 <= 1000+100 and $3 >= 1000-100){ ### excise region cover 100 bp of breakpoint
               print "SV Contig $contig: $6,$7\n";
               #$indel=1; ### indicate SV exists
               my $seqlen=length $refseq->{$contig};
               my $start=$6 < $7 ? $6 : $7;
               my $len  =abs($7-$6+1);
               $sequence=substr($refseq->{$contig},$start,$len);
               my $gap;
               while($sequence=~/(N+)/g){
                  $gap+=length $1;
               }
               if ($gap > $len*0.3){
                  $indel=3 unless ($indel == 1 or $indel == 0); ### too much unknown sequence in SV sequence, set to not sure
               }elsif ($len < 10){
                  $indel=0 unless ($indel == 1); ### insertion sequence too small, probably not a SV
               }else{
                  $indel=1;
               }
            }elsif(($1 <= 1000-100 and $2 >= 1000+100) or ($3 <= 1000-100 and $4 >= 1000+100)){
               $indel=0 unless ($indel == 1); ### indicate no SV exists
            }else{
               $indel=3 unless ($indel == 1 or $indel == 0); ### indicate not sure
            }
         }
      }
   }
}
close IN;
$/="\n";
return ($indel,$sequence);
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


