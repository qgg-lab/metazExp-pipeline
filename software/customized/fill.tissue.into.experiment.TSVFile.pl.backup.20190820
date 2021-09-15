#!/usr/bin/perl
use strict;
use DBI;
use Getopt::Long;
use DBI;
if($#ARGV<0){
        print "\nperl $0 \\\n" .
		"--inputExperimentTsvFile origin.experiment.tsv\\\n" .
		"--inputExperimentPrjIdFile prj.infor.tsv \\\n" .
		"--inputExperimentTissueFile manual.annot.with.tissue.tsv \\\n" .
		"--taxonId 9031 \\\n" .
		"--outputExperimentTsvFile experiment.with.tissue.and.prjId.tsv \n\n";
        exit;
}

my ($inputExperimentTsvFile, $inputExperimentPrjIdFile, $inputExperimentTissueFile, $taxonId, $outputExperimentTsvFile);
GetOptions(
	'inputExperimentTsvFile=s'=>\$inputExperimentTsvFile,
	'inputExperimentPrjIdFile=s'=>\$inputExperimentPrjIdFile,
	'inputExperimentTissueFile=s'=>\$inputExperimentTissueFile,
	'taxonId=s'=>\$taxonId,
	'outputExperimentTsvFile=s'=>\$outputExperimentTsvFile,
);

# read tissue information
my (%expInfo);
my ($expId, $tissue, $line);
open FF, "<$inputExperimentTissueFile";
#DRX001562       embryos
while($line = <FF>){
	chomp($line);
	($expId, $tissue) = split(/\t/, $line);
	${$expInfo{$expId}}{"tissue"} = $tissue;
}
close FF;

my ($prjId);
open FF, "<$inputExperimentPrjIdFile";
#DRX001555       DRP000595
while($line = <FF>){
	chomp($line);
	($expId, $prjId) = split(/\t/, $line);
	${$expInfo{$expId}}{"prjId"} = $prjId;
}
close FF;

my (@field);
my ($fieldNameString, $valueString);
open FF, "<$inputExperimentTsvFile";
# alignPercent, mappedSpots, experimentId, species, libraryType, libraryLayout, readLen, phredScore, totalSpots, runIdList, runNum, studyId, jcecTotalAsNum, jcecA5ssPercentage, jcecA3ssPercentage, jcecSePercentage, jcecRiPercentage, jcecMxePercentage, jcTotalAsNum, jcA5ssPercentage, jcA3ssPercentage, jcSePercentage, jcRiPercentage, jcMxePercentage_____94.41, 66.05, "SRX1036607", "Gallus gallus", "UN", "PAIRED", "101", "33", 69.97, "SRR2037196", 1, "SRP058621", 68536, 11.27, 13.67, 42.43, 28.47, 4.14, 68007, 11.25, 13.71, 42.47, 28.41, 4.13
open WW, ">$outputExperimentTsvFile";
while($line=<FF>){
	chomp($line);
	($fieldNameString, $valueString) = split(/_____/, $line);

	@field = split(/, /, $valueString);
	$expId = substr($field[2], 1, length($field[2]) - 2);

	print WW $fieldNameString . ", " . "tissue";
	print WW "_____";
	print WW $valueString . ", \"" . ${$expInfo{$expId}}{"tissue"} . "\"\n";
}
close FF;
close WW;
