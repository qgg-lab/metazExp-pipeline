#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--blastRlt blast.pairwise.txt \\\n" .
                "--outputDir  ./tmp \\\n" .
		"--searchId	11111  \\\n" .
                "--hitIsoformIdListFile \n\n\n";
	exit;
}

my ($blastRlt, $searchId, $outputDir, $hitIsoformIdListFile, $prefix);
GetOptions(
        'blastRlt=s'=>\$blastRlt,
        'outputDir=s'=>\$outputDir,
	'searchId=s'=>\$searchId,
        'hitIsoformIdListFile=s'=>\$hitIsoformIdListFile,
);

my ($blastText, $line, @hitIsoform, @hitIsoformIdList, $hitText, $isoformId);
open FF, "<$blastRlt";
while($line=<FF>){
	$blastText.=$line;
}
close FF;

$blastText = substr($blastText, 0, index($blastText, "Lambda"));
@hitIsoform = split(/>/, $blastText);
shift(@hitIsoform);

open WID, ">" . "$hitIsoformIdListFile";
foreach $hitText(@hitIsoform){
	if($hitText=~/^(.*?) .*/){
		$isoformId = $1;
		print WID $isoformId . "\n";
		open WW, ">" . $outputDir . "/blast.pairwise." . $isoformId . "." . $searchId;
		print WW ">" . $hitText;
		close WW;
	}
}
close WID;
