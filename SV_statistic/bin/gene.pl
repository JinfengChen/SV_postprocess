
my %data;
my $reflen=getfastalen("./HEG4_anno/HEG4.Deletion.final.fa");
my $cdslen=getfastalen("/rhome/cjinfeng/HEG4_cjinfeng/seqlib/GFF/MSU7/MSU7.all.cds");
my $cdsanno=readanno("/rhome/cjinfeng/HEG4_cjinfeng/seqlib/GFF/MSU7.gene.anno");
sumgene("./HEG4_anno/HEG4.Deletion.final.fa.genewise.gff",$reflen,$cdslen,$cdsanno,\%data);

#########
#1614   16.2  0.3  0.3  Chr10:10264822-10265112      1   291     (0) C SPMLIKE             DNA/En-Spm          (1911)   9078    8788    1
sub sumgene
{
my ($gff,$reflen,$cdslen,$cdsanno,$data)=@_;
my $refgff=parseGFFlen($gff);

print "SV\tSV_length\tTE_length\tSingle_Event\tElement:Rate\tElement:Family\n";
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
}

}


sub parseGFFlen
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;
my $id;    ##ID for element
my $record;##all line for this record
open IN, "$gff" or die "$!";
while (<IN>){
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
close IN;
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
    $hash{$head}=length $seq;
}
$/="\n";
return \%hash;
}




sub readanno
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}

