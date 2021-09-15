#!/usr/bin/perl
use strict;
use Getopt::Long;
use Bio::Seq;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--origOrfInCdna \\\n" .
                "--pepFasta \\\n" .
                "--cDNAFasta \\\n" .
                "--correctOrfInCda \n";
	exit;
}

my ($origOrfInCdna, $pepFasta, $cDNAFasta, $correctOrfInCda);

GetOptions(
        'origOrfInCdna=s'=>\$origOrfInCdna,
        'pepFasta=s'=>\$pepFasta,
        'cDNAFasta=s'=>\$cDNAFasta,
        'correctOrfInCda=s'=>\$correctOrfInCda,
);

# 将pep读入hash中
my (%pep,  $line, @seqId, $seqId);
open FF, "<$pepFasta";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$seqId = $1;
		@seqId = split(/ /, $seqId);
		$seqId = $seqId[0];
	}else{
		$pep{$seqId}.=$line;
	}
}
close FF;
my @pepId=keys(%pep);
print "finish load pep fasta:" . $#pepId . "\n";;

# 将cDNA读入hash中
my (%cDNAseq,  $line, @seqId, $seqId);
open FF, "<$cDNAFasta";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$seqId = $1;
		@seqId = split(/ /, $seqId);
		$seqId = $seqId[0];
	}else{
		$cDNAseq{$seqId}.=$line;
	}
}
close FF;
my @cDNAId = keys(%cDNAseq);
print "finish load cDNA fasta:" . $#cDNAId . "\n";

# 将orf读入，翻译orf，然后比较蛋白序列是否和hash中的蛋白序列一致
open FF, "<$origOrfInCdna";
open WW, ">$correctOrfInCda";
my ($orfBegin, $orfEnd, $orfSeq, @field);
# AT3G53400.1     420     1817    420     422     1818    1820
# AT5G45110.2     515     2062    515     517     2063    2065 
while($line=<FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);
	if(($field[2] - $field[1] + 1)/3 == length($pep{$field[0]})){
		print WW $line . "\n";
	}else{
		for(my $i=0; $i<=2; $i++){
			$orfBegin = $field[1] + $i;
			$orfEnd = length($pep{$field[0]})*3 + $orfBegin - 1;
			$orfSeq = uc(substr($cDNAseq{$field[0]}, $orfBegin-1, length($pep{$field[0]}) * 3));
			if(&TranslateDNASeq($orfSeq) eq $pep{$field[0]}){
				print WW join("\t", $field[0], $orfBegin, $orfEnd, $field[3], $field[4], $field[5], $field[6]) . "\n";
				last;
			}
		}
	}
}
close FF;
close WW;

sub TranslateDNASeq{
	use Bio::Seq;
	my ($dna)=@_;
	my $seqobj=Bio::Seq->new(-seq =>$dna, -alphabet =>'dna');
	return $seqobj->translate()->seq();
}


