#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;

my ($origAnnoName, $taxonId, $speciesName, $speciesAbbr, $genomeAssemblyVersion, $genomeAssemblyAddress, $genomeAnnoVersion, $genomeAnnoGtfAddress, $improveGtfMinMappedVolum, $improveGtfMinAligPer, $improveGtfExptUnitVolume, $improveGtfExptUnitTrsptMinCov, $improveGtfExptUnitExonMinCov, $improveGtfExptNumForImproveGtf, $improveGtfExptListFile, $identifyAsMinMappedVolum, $identifyAsMinAlignPercentage, $identifyAsExptNum, $identifyAsExptListFile, $origGtfGeneNum, $origGtfTrsptNum, $origAllAsNum, $origA3ssNum, $origA5ssNum, $origMxeNum, $origRiNum, $origSeNum, $improvedGtfTrsptNum, $improvedAllAsNum, $improvedA3ssNum, $improvedA5ssNum, $improvedMxeNum, $improvedRiNum, $improvedSeNum, $finalAllAsNum, $finalA3ssNum, $finalA5ssNum, $finalMxeNum, $finalRiNum, $finalSeNum, $origInnerConserAllAsNum, $origInnerConserA3ssNum, $origInnerConserA5ssNum, $origInnerConserMxeNum, $origInnerConserRiNum, $origInnerConserSeNum, $origOrthConserAllAsNum, $origOrthConserA3ssNum, $origOrthConserA5ssNum, $origOrthConserMxeNum, $origOrthConserRiNum, $origOrthConserSeNum, $improvedInnerConserAllAsNum, $improvedInnerConserA3ssNum, $improvedInnerConserA5ssNum, $improvedInnerConserMxeNum, $improvedInnerConserRiNum, $improvedInnerConserSeNum, $improvedOrthConserAllAsNum, $improvedOrthConserA3ssNum, $improvedOrthConserA5ssNum, $improvedOrthConserMxeNum, $improvedOrthConserRiNum, $improvedOrthConserSeNum, $finalInnerConserAllAsNum, $finalInnerConservedA3ssNum, $finalInnerConservedA5ssNum, $finalInnerConservedMxeNum, $finalInnerConservedRiNum, $finalInnerConservedSeNum, $finalOrthConserAllAsNum, $finalOrthConserA3ssNum, $finalOrthConserA5ssNum, $finalOrthConserMxeNum, $finalOrthConserRiNum, $finalOrthConserSeNum, $speciesDesptText, $annoDesptText, $identifyAsDesptText, $origExonNum, $origSpliceNum, $improvedExonNum, $improvedSpliceNum, $outputSpeciesTsv, $collectedAllExptListFile, $collectedAllExptNum);

GetOptions(
	'speciesDesptText=s'=>\$speciesDesptText,
	'annoDesptText=s'=>\$annoDesptText,
	'identifyAsDesptText=s'=>\$identifyAsDesptText,
	'origAnnoName=s'=>\$origAnnoName,
	'taxonId=s'=>\$taxonId,
	'speciesName=s'=>\$speciesName,
	'speciesAbbr=s'=>\$speciesAbbr,
	'genomeAssemblyVersion=s'=>\$genomeAssemblyVersion,
	'genomeAssemblyAddress=s'=>\$genomeAssemblyAddress,
	'genomeAnnoVersion=s'=>\$genomeAnnoVersion,
	'genomeAnnoGtfAddress=s'=>\$genomeAnnoGtfAddress,
	'improveGtfMinMappedVolum=s'=>\$improveGtfMinMappedVolum,
	'improveGtfMinAligPer=s'=>\$improveGtfMinAligPer,
	'improveGtfExptUnitVolume=s'=>\$improveGtfExptUnitVolume,
	'improveGtfExptUnitTrsptMinCov=s'=>\$improveGtfExptUnitTrsptMinCov,
	'improveGtfExptUnitExonMinCov=s'=>\$improveGtfExptUnitExonMinCov,
	'improveGtfExptNumForImproveGtf=s'=>\$improveGtfExptNumForImproveGtf,
	'improveGtfExptListFile=s'=>\$improveGtfExptListFile,
	'identifyAsMinMappedVolum=s'=>\$identifyAsMinMappedVolum,
	'identifyAsMinAlignPercentage=s'=>\$identifyAsMinAlignPercentage,
	'identifyAsExptNum=s'=>\$identifyAsExptNum,
	'collectedAllExptListFile=s'=>\$collectedAllExptListFile,
	'collectedAllExptNum=s'=>\$collectedAllExptNum,
	'identifyAsExptListFile=s'=>\$identifyAsExptListFile,
	'origGtfGeneNum=s'=>\$origGtfGeneNum,
	'origGtfTrsptNum=s'=>\$origGtfTrsptNum,
	'origAllAsNum=s'=>\$origAllAsNum,
	'origA3ssNum=s'=>\$origA3ssNum,
	'origA5ssNum=s'=>\$origA5ssNum,
	'origMxeNum=s'=>\$origMxeNum,
	'origRiNum=s'=>\$origRiNum,
	'origSeNum=s'=>\$origSeNum,
	'improvedGtfTrsptNum=s'=>\$improvedGtfTrsptNum,
	'improvedAllAsNum=s'=>\$improvedAllAsNum,
	'improvedA3ssNum=s'=>\$improvedA3ssNum,
	'improvedA5ssNum=s'=>\$improvedA5ssNum,
	'improvedMxeNum=s'=>\$improvedMxeNum,
	'improvedRiNum=s'=>\$improvedRiNum,
	'improvedSeNum=s'=>\$improvedSeNum,
	'finalAllAsNum=s'=>\$finalAllAsNum,
	'finalA3ssNum=s'=>\$finalA3ssNum,
	'finalA5ssNum=s'=>\$finalA5ssNum,
	'finalMxeNum=s'=>\$finalMxeNum,
	'finalRiNum=s'=>\$finalRiNum,
	'finalSeNum=s'=>\$finalSeNum,
	'origInnerConserAllAsNum=s'=>\$origInnerConserAllAsNum,
	'origInnerConserA3ssNum=s'=>\$origInnerConserA3ssNum,
	'origInnerConserA5ssNum=s'=>\$origInnerConserA5ssNum,
	'origInnerConserMxeNum=s'=>\$origInnerConserMxeNum,
	'origInnerConserRiNum=s'=>\$origInnerConserRiNum,
	'origInnerConserSeNum=s'=>\$origInnerConserSeNum,
	'origOrthConserAllAsNum=s'=>\$origOrthConserAllAsNum,
	'origOrthConserA3ssNum=s'=>\$origOrthConserA3ssNum,
	'origOrthConserA5ssNum=s'=>\$origOrthConserA5ssNum,
	'origOrthConserMxeNum=s'=>\$origOrthConserMxeNum,
	'origOrthConserRiNum=s'=>\$origOrthConserRiNum,
	'origOrthConserSeNum=s'=>\$origOrthConserSeNum,
	'improvedInnerConserAllAsNum=s'=>\$improvedInnerConserAllAsNum,
	'improvedInnerConserA3ssNum=s'=>\$improvedInnerConserA3ssNum,
	'improvedInnerConserA5ssNum=s'=>\$improvedInnerConserA5ssNum,
	'improvedInnerConserMxeNum=s'=>\$improvedInnerConserMxeNum,
	'improvedInnerConserRiNum=s'=>\$improvedInnerConserRiNum,
	'improvedInnerConserSeNum=s'=>\$improvedInnerConserSeNum,
	'improvedOrthConserAllAsNum=s'=>\$improvedOrthConserAllAsNum,
	'improvedOrthConserA3ssNum=s'=>\$improvedOrthConserA3ssNum,
	'improvedOrthConserA5ssNum=s'=>\$improvedOrthConserA5ssNum,
	'improvedOrthConserMxeNum=s'=>\$improvedOrthConserMxeNum,
	'improvedOrthConserRiNum=s'=>\$improvedOrthConserRiNum,
	'improvedOrthConserSeNum=s'=>\$improvedOrthConserSeNum,
	'finalInnerConserAllAsNum=s'=>\$finalInnerConserAllAsNum,
	'finalInnerConservedA3ssNum=s'=>\$finalInnerConservedA3ssNum,
	'finalInnerConservedA5ssNum=s'=>\$finalInnerConservedA5ssNum,
	'finalInnerConservedMxeNum=s'=>\$finalInnerConservedMxeNum,
	'finalInnerConservedRiNum=s'=>\$finalInnerConservedRiNum,
	'finalInnerConservedSeNum=s'=>\$finalInnerConservedSeNum,
	'finalOrthConserAllAsNum=s'=>\$finalOrthConserAllAsNum,
	'finalOrthConserA3ssNum=s'=>\$finalOrthConserA3ssNum,
	'finalOrthConserA5ssNum=s'=>\$finalOrthConserA5ssNum,
	'finalOrthConserMxeNum=s'=>\$finalOrthConserMxeNum,
	'finalOrthConserRiNum=s'=>\$finalOrthConserRiNum,
	'finalOrthConserSeNum=s'=>\$finalOrthConserSeNum,
	'origExonNum=s'=>\$origExonNum,
	'origSpliceNum=s'=>\$origSpliceNum,
	'improvedExonNum=s'=>\$improvedExonNum,
	'improvedSpliceNum=s'=>\$improvedSpliceNum,
	'outputSpeciesTsv=s'=>\$outputSpeciesTsv
);

# ???????????????gtf?????????experiment??????
my ($exptId, $improvedGtfExptList, $identifyAsExptList);
open FF, "<$improveGtfExptListFile";
while($exptId=<FF>){
	chomp($exptId);
	$improvedGtfExptList .= $exptId . ",";
}
close FF;
if($improvedGtfExptList ne ""){
	$improvedGtfExptList = substr($improvedGtfExptList, 0, length($improvedGtfExptList) - 1);
}else{
	$improvedGtfExptList = "none";
}

open FF, "<$identifyAsExptListFile";
while($exptId=<FF>){
	chomp($exptId);
	$identifyAsExptList .= $exptId . ",";
}
close FF;
if($identifyAsExptList ne ""){
	$identifyAsExptList = substr($identifyAsExptList, 0, length($identifyAsExptList) - 1);
}else{
	$identifyAsExptList = "none";
}
my $collectedAllExptList="";
open FF, "<$collectedAllExptListFile";
while($exptId=<FF>){
	chomp($exptId);
	$collectedAllExptList .= $exptId . ",";
}
close FF;
$collectedAllExptList = substr($collectedAllExptList, 0, length($collectedAllExptList) - 1);

my ($fieldNameString, $fieldValueString);
$fieldNameString = join("###", "taxonId", "speciesName", "speciesAbbr", "genomeAssemblyVersion", "genomeAssemblyAddress", "genomeAnnoVersion", "genomeAnnoGtfAddress", "improveGtfMinMappedVolum", "improveGtfMinAligPer", "improveGtfExptUnitVolume", "improveGtfExptUnitTrsptMinCov", "improveGtfExptUnitExonMinCov", "improveGtfExptNumForImproveGtf", "identifyAsMinMappedVolum", "identifyAsMinAlignPercentage", "identifyAsExptNum", "collectedAllExptNum", "origGtfGeneNum", "origGtfTrsptNum", "origAllAsNum", "origA3ssNum", "origA5ssNum", "origMxeNum", "origRiNum", "origSeNum", "improvedGtfTrsptNum", "improvedAllAsNum", "improvedA3ssNum", "improvedA5ssNum", "improvedMxeNum", "improvedRiNum", "improvedSeNum", "finalAllAsNum", "finalA3ssNum", "finalA5ssNum", "finalMxeNum", "finalRiNum", "finalSeNum", "origInnerConserAllAsNum", "origInnerConserA3ssNum", "origInnerConserA5ssNum", "origInnerConserMxeNum", "origInnerConserRiNum", "origInnerConserSeNum", "origOrthConserAllAsNum", "origOrthConserA3ssNum", "origOrthConserA5ssNum", "origOrthConserMxeNum", "origOrthConserRiNum", "origOrthConserSeNum", "improvedInnerConserAllAsNum", "improvedInnerConserA3ssNum", "improvedInnerConserA5ssNum", "improvedInnerConserMxeNum", "improvedInnerConserRiNum", "improvedInnerConserSeNum", "improvedOrthConserAllAsNum", "improvedOrthConserA3ssNum", "improvedOrthConserA5ssNum", "improvedOrthConserMxeNum", "improvedOrthConserRiNum", "improvedOrthConserSeNum", "finalInnerConserAllAsNum", "finalInnerConservedA3ssNum", "finalInnerConservedA5ssNum", "finalInnerConservedMxeNum", "finalInnerConservedRiNum", "finalInnerConservedSeNum", "finalOrthConserAllAsNum", "finalOrthConserA3ssNum", "finalOrthConserA5ssNum", "finalOrthConserMxeNum", "finalOrthConserRiNum", "finalOrthConserSeNum", "improvedGtfExptList", "identifyAsExptList", "collectedAllExptList", "origAnnoName", "speciesDesptText", "annoDesptText", "identifyAsDesptText", "origExonNum", "origSpliceNum", "improvedExonNum", "improvedSpliceNum");

$speciesDesptText =~s/ /__/g;
$annoDesptText =~s/ /__/g;
$identifyAsDesptText=~s/ /__/g;

$fieldValueString = join("###", $taxonId, $speciesName, $speciesAbbr, $genomeAssemblyVersion, $genomeAssemblyAddress, $genomeAnnoVersion, $genomeAnnoGtfAddress, $improveGtfMinMappedVolum, $improveGtfMinAligPer, $improveGtfExptUnitVolume, $improveGtfExptUnitTrsptMinCov, $improveGtfExptUnitExonMinCov, $improveGtfExptNumForImproveGtf, $identifyAsMinMappedVolum, $identifyAsMinAlignPercentage, $identifyAsExptNum, $collectedAllExptNum, $origGtfGeneNum, $origGtfTrsptNum, $origAllAsNum, $origA3ssNum, $origA5ssNum, $origMxeNum, $origRiNum, $origSeNum, $improvedGtfTrsptNum, $improvedAllAsNum, $improvedA3ssNum, $improvedA5ssNum, $improvedMxeNum, $improvedRiNum, $improvedSeNum, $finalAllAsNum, $finalA3ssNum, $finalA5ssNum, $finalMxeNum, $finalRiNum, $finalSeNum, $origInnerConserAllAsNum, $origInnerConserA3ssNum, $origInnerConserA5ssNum, $origInnerConserMxeNum, $origInnerConserRiNum, $origInnerConserSeNum, $origOrthConserAllAsNum, $origOrthConserA3ssNum, $origOrthConserA5ssNum, $origOrthConserMxeNum, $origOrthConserRiNum, $origOrthConserSeNum, $improvedInnerConserAllAsNum, $improvedInnerConserA3ssNum, $improvedInnerConserA5ssNum, $improvedInnerConserMxeNum, $improvedInnerConserRiNum, $improvedInnerConserSeNum, $improvedOrthConserAllAsNum, $improvedOrthConserA3ssNum, $improvedOrthConserA5ssNum, $improvedOrthConserMxeNum, $improvedOrthConserRiNum, $improvedOrthConserSeNum, $finalInnerConserAllAsNum, $finalInnerConservedA3ssNum, $finalInnerConservedA5ssNum, $finalInnerConservedMxeNum, $finalInnerConservedRiNum, $finalInnerConservedSeNum, $finalOrthConserAllAsNum, $finalOrthConserA3ssNum, $finalOrthConserA5ssNum, $finalOrthConserMxeNum, $finalOrthConserRiNum, $finalOrthConserSeNum, $improvedGtfExptList, $identifyAsExptList, $collectedAllExptList, $origAnnoName, $speciesDesptText, $annoDesptText, $identifyAsDesptText, $origExonNum, $origSpliceNum, $improvedExonNum, $improvedSpliceNum);

open WW, ">$outputSpeciesTsv";
print WW $fieldNameString . "___" . $fieldValueString;
close WW;
