#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--inputExpTsv original.livestock.experiment.tsv\\\n" .
                "--inputExpRawSraInfo total.exp.raw.sra.info.txt \\\n" .
                "--outputExpTsv livestock.experiment.tsv \n";
	exit;
}

my ($inputExpTsv, $inputExpRawSraInfo, $outputExpTsv);

GetOptions(
        'inputExpTsv=s'=>\$inputExpTsv,
        'inputExpRawSraInfo=s'=>\$inputExpRawSraInfo,
        'outputExpTsv=s'=>\$outputExpTsv,
);


my ($line, @field, $field, $expId, $studyTitle, $studyDescription, $experimentTitle, $experimentDescription);
my (%exptInfo);
open FF, "<$inputExpRawSraInfo";
while($line = <FF>){
	chomp($line);
	@field = ();
	@field = split(/\|___\|/, $line);

	$expId = $field[18];
	$studyTitle = $field[58];
	$studyTitle=~s/"/\\"/g;
	$studyTitle=~s/'/\\'/g;
	$studyDescription = $field[62];
	$studyDescription =~s/"/\\"/g;
	$studyDescription =~s/'/\\'/g;
	$experimentTitle = $field[19];
	$experimentTitle =~s/"/\\"/g;
	$experimentTitle =~s/'/\\'/g;
	${$exptInfo{$expId}}{"studyTitle"} = $studyTitle;
	${$exptInfo{$expId}}{"studyDescription"} = $studyDescription;
	${$exptInfo{$expId}}{"exptTitle"} = $experimentTitle;
}
close FF;

my ($fieldNameString, $valueString);
open FF, "<$inputExpTsv";
# alignPercent, mappedSpots, experimentId, species, libraryType, libraryLayout, readLen, phredScore, totalSpots, runIdList, runNum, studyId, jcecTotalAsNum, jcecA5ssPercentage, jcecA3ssPercentage, jcecSePercentage, jcecRiPercentage, jcecMxePercentage, jcTotalAsNum, jcA5ssPercentage, jcA3ssPercentage, jcSePercentage, jcRiPercentage, jcMxePercentage, tissue_____94.41, 66.05, "SRX1036607", "Gallus gallus", "UN", "PAIRED", "101", "33", 69.97, "SRR2037196", 1, "SRP058621", 68536, 11.27, 13.67, 42.43, 28.47, 4.14, 68007, 11.25, 13.71, 42.47, 28.41, 4.13, "embryos" 
open WW, ">$outputExpTsv";
while($line=<FF>){
        chomp($line);
        ($fieldNameString, $valueString) = split(/_____/, $line);

        @field = split(/, /, $valueString);
        $expId = substr($field[2], 1, length($field[2]) - 2);

        print WW $fieldNameString . ", studyTitle, studyDesc, exptTitle";
        print WW "_____";
        print WW $valueString . ", \"" . ${$exptInfo{$expId}}{"studyTitle"} . "\", \"" . ${$exptInfo{$expId}}{"studyDescription"} . "\", \"" . ${$exptInfo{$expId}}{"exptTitle"} . "\"\n";
}
close FF;
close WW;

