#!/usr/bin/perl
use strict;
use Getopt::Long;

if($#ARGV<0){
	print "\nperl $0 \\\n" .
		"\t\t--mappedExonFileList 9031.bed,9796.bed,9823.bed,9913.bed,9940.bed \\\n" .
		"\t\t--taxonIdList 9031,9796,9823,9913,9940 \\\n" .
		"\t\t--outputSingleCopyExon orthExon.tsv \\\n" .
		"\t\t--outputTotalExon totalExon.tsv \n";
	exit;
}

my ($mappedExonFileList, $outputSingleCopyExon, $taxonIdList, $outputTotalExon);

GetOptions(
        'mappedExonFileList=s'=>\$mappedExonFileList,
	'taxonIdList=s'=>\$taxonIdList,
        'outputSingleCopyExon=s'=>\$outputSingleCopyExon,
	'outputTotalExon=s'=>\$outputTotalExon,
);

my (@taxonId, $taxonId);
my (@mappedExonFile, $mappedExonFile, $line, @fields, $exonId, $chain, $chr, $start, $stop, $orthId, $score);
my (%orthIdToExonId, %orthIdToExonNum);
my (@orthId, $orthId, $flag);
my ($orthNum, $orthTag);
my ($totalNum);

@taxonId = split(/,/, $taxonIdList);
@mappedExonFile = split(/,/, $mappedExonFileList);

for(my $i=0; $i<=$#taxonId; $i++){
	$mappedExonFile = $mappedExonFile[$i];
	$taxonId = $taxonId[$i];
	open FF, "<$mappedExonFile";
	while($line=<FF>){
		chomp($line);
		($chr, $start, $stop, $exonId, $score, $chain) = ("", "", "", "", "", "");
		($chr, $start, $stop, $exonId, $score, $chain) = split(/\t/, $line);
		$orthId = join("#", $chr, $start, $stop, $chain);

		${$orthIdToExonId{$orthId}}{$taxonId}.=$exonId . "#";
		${$orthIdToExonNum{$orthId}}{$taxonId}++;
	}
	close FF;
}


# 输出
open WT, ">$outputTotalExon";
open WW, ">$outputSingleCopyExon";
print WW "orthId\torthExonCoord";
print WT "orthId";
foreach $taxonId(@taxonId){
	print WW "\t" . $taxonId;
	print WT "\t" . $taxonId;
}
print WW "\n";
print WT "\n";

@orthId = keys(%orthIdToExonId);

foreach $orthId(@orthId){
	$flag = 1;
	$orthNum++;
	$orthTag = "ORTH" . sprintf("%010d", $orthNum);

	print WT $orthTag;
	foreach $taxonId(@taxonId){
		$flag = 0 if(${$orthIdToExonNum{$orthId}}{$taxonId}!=1);
		print WT "\t" . ${$orthIdToExonNum{$orthId}}{$taxonId};
	}
	print WT "\n";

	next if($flag ==0);
	
	print WW $orthTag . "\t" . $orthId;
	foreach $taxonId(@taxonId){
		${$orthIdToExonId{$orthId}}{$taxonId} = substr(${$orthIdToExonId{$orthId}}{$taxonId}, 0, length(${$orthIdToExonId{$orthId}}{$taxonId}) -1);
		print WW "\t" . ${$orthIdToExonId{$orthId}}{$taxonId};
	}

	print WW "\n";
}
close WW;
close WT;
