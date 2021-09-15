#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--inputOriginExonCoord exon.origin.coord.tsv\\\n" .
                "--inputOrthExonMatrix final.orthExon.matrix.tsv \\\n" .
                "--inputTaxonId 9031 \\\n" .
		"--outputOriginExonCoord orth.exon.coord.tsv \\\n";
	exit;
}

my ($inputOriginExonCoord, $inputOrthExonMatrix, $inputTaxonId, $outputOriginExonCoord);
GetOptions(
        'inputOriginExonCoord=s'=>\$inputOriginExonCoord,
        'inputOrthExonMatrix=s'=>\$inputOrthExonMatrix,
        'inputTaxonId=s'=>\$inputTaxonId,
        'outputOriginExonCoord=s'=>\$outputOriginExonCoord,
);

my (%orthId);
my (@fields, $line, @taxonId, $taxonId, $orthId, $i, $taxonPosition);

my (%originExonCoord);
open FF, "<$inputOriginExonCoord";
#chr1    5273    5524    Exon1   0       -
while($line=<FF>){
	chomp($line);
	@fields = split(/\t/, $line);
	$originExonCoord{$fields[3]}=$line;
}
close FF;

open WW, ">$outputOriginExonCoord";
open FF, "<$inputOrthExonMatrix";
#orthId  9031    9796    9823    9913    9940
#ORTH0000000002  Exon116411      Exon225128      Exon118201      Exon213018      Exon185323
$line = <FF>;
chomp($line);
@taxonId = split(/\t/, $line);
shift(@taxonId);
$taxonPosition = -1;
for($i=0; $i<=$#taxonId; $i++){
	if($taxonId[$i] eq $inputTaxonId){
		$taxonPosition = $i;
		last;
	}
}

while($line=<FF>){
	chomp($line);
	@fields = ();
	@fields = split(/\t/, $line);
	$orthId = shift(@fields);
	print WW join("\t", $originExonCoord{$fields[$taxonPosition]}) . "\t" . $orthId . "\n";
}
close FF;
close WW;
