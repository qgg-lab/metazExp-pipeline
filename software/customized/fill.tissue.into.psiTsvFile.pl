#!/usr/bin/perl
use strict;
use DBI;
use Getopt::Long;
use DBI;
if($#ARGV<0){
        print "\nperl $0 \\\n" .
		"--inputPsiTsvFile psi.tsv \\\n" .
		"--inputExpInformation experiment.with.tissue.and.prjId.tsv \\\n" .
		"--outputPsiTsvFile new.psi.tsv \n\n";
        exit;
}

my ($inputExpInformation, $inputPsiTsvFile, $outputPsiTsvFile);
GetOptions(
	'inputExpInformation=s'=>\$inputExpInformation,
	'inputPsiTsvFile=s'=>\$inputPsiTsvFile,
	'outputPsiTsvFile=s'=>\$outputPsiTsvFile,
);

my (@fieldTitle, $fieldTitle, @field, $field, $line, @tmp);

my ($fieldNameString, $valueString, $expId, $tissue, %expIdToTissue);
# 读取experiment信息，建立expId -> tissue 关系
open FF, "<$inputExpInformation";
while($line = <FF>){
	chomp($line);
	($fieldNameString, $valueString) = split(/_____/, $line);

	@field = split(/, /, $valueString);
	$expId = substr($field[2], 1, length($field[2]) - 2);
	$tissue = substr($field[$#field], 1, length($field[$#field]) - 2);
	$expIdToTissue{$expId} = $tissue;
}
close FF;

# 读取psi,填充 experiment -> tissue

open FF, "<$inputPsiTsvFile";
open WW, ">$outputPsiTsvFile";
while($line = <FF>){
	chomp($line);
	($fieldNameString, $valueString) = split(/_____/, $line);

	@field = split(/, /, $valueString);

	$expId = substr($field[2], 1, length($field[2]) - 2);
	print WW $fieldNameString . ", " . "tissue" . "_____";
	print WW $valueString . ", \"" . $expIdToTissue{$expId} . "\"\n";
}
close FF;


