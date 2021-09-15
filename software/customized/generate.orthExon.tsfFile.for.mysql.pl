#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--exonFileList 9031.origin.exon.bed,9796.origin.exon.bed,9823.origin.exon.bed,9913.origin.exon.bed,/9940.origin.exon.bed \\\n" .
		"--taxonList 9031,9796,9823,9913,9940 \n" .
		"--orthExonMatrix final.orthExon.matrix.tsv \n" .
		"--outputOrthExonTsv orthExon.tsv\\\n";
	exit;
}

my ($exonFileList, $taxonList, $orthExonMatrix, $outputOrthExonTsv);
my (@exonFile, $exonFile, @taxon, $taxon);
GetOptions(
        'exonFileList=s'=>\$exonFileList,
	'taxonList=s'=>\$taxonList,
	'orthExonMatrix=s'=>\$orthExonMatrix, 
	'outputOrthExonTsv=s'=>\$outputOrthExonTsv,
);

my (@fieldTitle, $fieldTitle, @field, $field, $line, $orthExonId, %exonToOrthExonId);

open FF, "<$orthExonMatrix";
$line = <FF>;
chomp($line);
@fieldTitle = split(/\t/, $line);
shift(@fieldTitle);
# orthId  9031    9796    9823    9913    9940
# ORTH0000000002  Exon116411      Exon225128      Exon118201      Exon213018      Exon185323
while($line=<FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);
	$orthExonId = shift(@field);
	for(my $i=0; $i<=$#field; $i++){
		$exonToOrthExonId{$fieldTitle[$i] . "#" . $field[$i]} = $orthExonId;
	}
}
close FF;

my ($chr);

open WW, ">$outputOrthExonTsv";
@exonFile = split(/,/, $exonFileList);
@taxon = split(/,/, $taxonList);

for(my $i=0; $i<=$#taxon; $i++){
	$exonFile = $exonFile[$i];
	$taxon = $taxon[$i];
	open FF, "<$exonFile";
	#chr1    5272    5524    Exon1   0       -
	while($line=<FF>){
		chomp($line);
		@field = ();
		@field = split(/\t/, $line);
		$chr = $field[0];
		if($chr=~/chr0*(.*)/){
			$chr = $1;
		}
		$orthExonId = $taxon . "#" . $field[3];
		if(exists($exonToOrthExonId{$orthExonId})){
			print WW join(", ", "taxonId", "orthExonId", "chr", "strand", "start", "stop");
			print WW "_____";
			print WW "\"" . $taxon . "\", " . "\"" . $exonToOrthExonId{$orthExonId} . "\", " .  "\"" . $chr . "\", \"" . $field[5] .  "\", " . $field[1] . ", " . $field[2] . "\n";
		}
	}
	close FF;
}	

close WW;
