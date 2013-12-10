
my %data;
my $reflen=getfastalen("./HEG4_anno/HEG4.insertion.final.fa");
sumrepeat("./HEG4_anno/HEG4.insertion.final.fa.out",$reflen,\%data);


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
     }
}

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

