#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--inputAsFile A5SS.catalog\\\n" .
		"--inputAsWithOrthExonFile A5SS.catalog.with.orthTag \\\n" .
		"--ASType A5SS \\\n" .
		"--outputCoordBedFile A5SS.coords.in.orthExon.bed \\\n";
	exit;
}

my ($inputAsFile, $inputAsWithOrthExonFile, $ASType, $outputCoordBedFile);
GetOptions(
	'inputAsFile=s'=>\$inputAsFile, 
	'inputAsWithOrthExonFile=s'=>\$inputAsWithOrthExonFile, 
	'ASType=s'=>\$ASType, 
	'outputCoordBedFile=s'=>\$outputCoordBedFile, 
);

my ($line, @fields, $field);

# read AS
my (@fieldTitle, $i,  %as, $orthId);
open FF, "<$inputAsFile";
#ASID GeneID geneSymbol chr strand longExonStart_0base longExonEnd shortES shortEE flankingES flankingEE
$line=<FF>;
chomp($line);
@fieldTitle = ();
@fieldTitle = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@fields = ();
	@fields = split(/\t/, $line);
	for($i=0; $i<=$#fields; $i++){
		${$as{$fields[0]}}{$fieldTitle[$i]}=$fields[$i];
	}		
}

my (%tmpAs);
open FF, "<$inputAsWithOrthExonFile";
$line=<FF>;
#ASID GeneID geneSymbol chr strand longExonStart_0base longExonEnd shortES shortEE flankingES flankingEE
chomp($line);
@fieldTitle = ();
@fieldTitle = split(/\t/, $line);
open WW, ">$outputCoordBedFile";
# GGALA5SS0000018079
# chr26    5272    5273    GGALA5SS0000018079:longExonStart_0base:ORTH0000196814   0       -
# chr26    7279    7280    GGALA5SS0000018079:longExonEnd:ORTH0000196814           0       -
# chr5     8971    8972    GGALA5SS0000009193:flankingES:ORTH0000196814            0	    -
# chr5     9931    9932    GGALA5SS0000009193:flankingEE:ORTH0000196814            0	    -
while($line=<FF>){
	#GGALA5SS0000018079 "ENSGALG00000000164" "MYBPH"    chr26 -      ORTH0000196814      ORTH0000196814 ORTH0000196814 ORTH0000196814 -          -
	#GGALA5SS0000009193 "ENSGALG00000000296" "PCNX1" chr5 - - - - - ORTH0000227401 ORTH0000227401
	chomp($line);
	@fields = ();
	@fields = split(/\t/, $line);
	for($i=5; $i<=$#fields; $i++){
		if($fields[$i]=~/ORTH/){
			if($fieldTitle[$i]=~/End/ or $fieldTitle[$i] =~/EE/){
				print WW join("\t", $fields[3], ${$as{$fields[0]}}{$fieldTitle[$i]}-1, ${$as{$fields[0]}}{$fieldTitle[$i]}, $fields[0] . ":" . $fields[$i] . ":" . $fieldTitle[$i], 0, $fields[4]) . "\n";
			}else{
				print WW join("\t", $fields[3], ${$as{$fields[0]}}{$fieldTitle[$i]}, ${$as{$fields[0]}}{$fieldTitle[$i]}+1, $fields[0] . ":" . $fields[$i] . ":" . $fieldTitle[$i], 0, $fields[4]) . "\n";
			}
		}
	}
}
close WW;
close FF;
