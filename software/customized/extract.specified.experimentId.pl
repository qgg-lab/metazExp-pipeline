#!/usr/bin/perl
use strict;
use Getopt::Long;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--curatedSampleTsv curated.experiment.tsv \\\n" .
                "--specifiedField  Assemble\\\n" .
		"--specifiedValue 1 \\\n" .
		"--outputExpIdFile outputExptId.list \n";
	exit;
}

my ($curatedSampleTsv, $specifiedField, $specifiedValue, $outputExpIdFile);

GetOptions(
        'curatedSampleTsv=s'=>\$curatedSampleTsv,
        'specifiedField=s'=>\$specifiedField,
        'specifiedValue=s'=>\$specifiedValue,
	'outputExpIdFile=s'=>\$outputExpIdFile,
);

# Tissue  SubTissue       TissueGroup     Development     Treatment       TreatmentGroup  Experiment      Study   Base    Layout  SpotsNum        ReadNum SpotLen ReadLen Gather  AS      Assemble        RunList Phenotype
my ($line, $i);
my (@fieldName, @field, $expIdPos, $specifiedFilePos);
#print join("\t", $specifiedField, $specifiedValue) . "\n";

open FF, "<$curatedSampleTsv";
open WW, ">$outputExpIdFile";
$line =<FF>;
chomp($line);
@fieldName = split(/\t/, $line);

for($i=0; $i<=$#fieldName; $i++){
	$expIdPos = $i if($fieldName[$i] eq "Experiment");
	$specifiedFilePos = $i if($fieldName[$i] eq $specifiedField);
}
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	print WW $field[$expIdPos] . "\n" if($field[$specifiedFilePos] eq $specifiedValue);
}

close FF;
close WW;
