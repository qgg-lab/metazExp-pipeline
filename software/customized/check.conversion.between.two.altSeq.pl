#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--asTsvMysql \\\n" .
		"--outputCheckRlt \n";
	exit;
}

my ($asTsvMysql, $outputCheckRlt);

GetOptions(
        'asTsvMysql=s'=>\$asTsvMysql,
        'outputCheckRlt=s'=>\$outputCheckRlt,
);

# 读取asTsvMysql
my (@nameField, @valueField, $nameFieldString, $valueFieldString);
my (%tmpHash, $i, $newExclSeq, $newInclSeq, $registedSeq, $preSubSeq, $tailSubSeq, $insertSubSeq);
my ($inclusionAltSeq, $exclusionAltSeq);
my ($editSiteInBothAltSeq, $cutSizeInInclAlt, $insertSizeInInclAlt, $cutSizeInExclAlt, $insertSizeInExclAlt);
open WW, ">$outputCheckRlt";
open FF, "<$asTsvMysql";
while(my $line=<FF>){
	chomp($line);
	($nameFieldString, $valueFieldString) = split(/___/, $line);
	@nameField = split(/, /, $nameFieldString);
	@valueField = split(/, /, $valueFieldString);
	for($i=0; $i<=$#nameField; $i++){
		$tmpHash{$nameField[$i]} = $valueField[$i];
	}
	
	($inclusionAltSeq, $exclusionAltSeq, $editSiteInBothAltSeq, $cutSizeInInclAlt, $insertSizeInInclAlt, $cutSizeInExclAlt, $insertSizeInExclAlt) = ($tmpHash{"inclusionAltSeq"}, $tmpHash{"exclusionAltSeq"}, $tmpHash{"editSiteInBothAltSeq"}, $tmpHash{"cutSizeInInclAlt"}, $tmpHash{"insertSizeInInclAlt"}, $tmpHash{"cutSizeInExclAlt"}, $tmpHash{"insertSizeInExclAlt"});

	# 检查从inclusionAltSeq到exclusionAltSeq
	$preSubSeq = substr($inclusionAltSeq, 0, $editSiteInBothAltSeq);
	$tailSubSeq = substr($inclusionAltSeq, $editSiteInBothAltSeq + $cutSizeInInclAlt);
	$insertSubSeq = substr($exclusionAltSeq, $editSiteInBothAltSeq, $cutSizeInExclAlt);
	if($tmpHash{"asType"} eq "A3SS" or $tmpHash{"asType"} eq "A5SS" or $tmpHash{"asType"} eq "RI" or $tmpHash{"asType"}  eq "SE"){
		$newExclSeq = $preSubSeq . $tailSubSeq;
	}else{
		$newExclSeq = $preSubSeq . $insertSubSeq . $tailSubSeq;
	}
	if($newExclSeq ne $exclusionAltSeq){
		print WW join("\t", $tmpHash{"asId"}, "error: inclusionAltSeq --> exclusionAltSeq") . "\n";
		print WW join("\t", "editSiteInBothAltSeq", "cutSizeInInclAlt", "insertSizeInInclAlt", "cutSizeInExclAlt", "insertSizeInExclAlt") . "\n";
		print WW join("\t", $editSiteInBothAltSeq, $cutSizeInInclAlt, $insertSizeInInclAlt, $cutSizeInExclAlt, $insertSizeInExclAlt) . "\n";
		print WW "inclusion:\n" . $inclusionAltSeq . "\n";
		print WW "exclusion:\n" . $exclusionAltSeq . "\n";
		print WW "newExclSeq:\n" . $newExclSeq . "\n";
	}

	# 检查从exclusionAltSeq到inclusionAltSeq
	$preSubSeq = substr($exclusionAltSeq, 0, $editSiteInBothAltSeq);
	$tailSubSeq = substr($exclusionAltSeq, $editSiteInBothAltSeq + $cutSizeInExclAlt);
	$insertSubSeq = substr($inclusionAltSeq, $editSiteInBothAltSeq, $cutSizeInInclAlt);
	$newInclSeq = $preSubSeq . $insertSubSeq . $tailSubSeq;
	
	if($newInclSeq ne $inclusionAltSeq){
		print WW join("\t", $tmpHash{"asId"}, "error: exclusionAltSeq --> inclusionAltSeq") . "\n";
		print WW join("\t", "editSiteInBothAltSeq", "cutSizeInInclAlt", "insertSizeInInclAlt", "cutSizeInExclAlt", "insertSizeInExclAlt") . "\n";
		print WW join("\t", $editSiteInBothAltSeq, $cutSizeInInclAlt, $insertSizeInInclAlt, $cutSizeInExclAlt, $insertSizeInExclAlt) . "\n";
		print WW "exclusion:\n" . $exclusionAltSeq . "\n";
		print WW "inclusion:\n" . $inclusionAltSeq . "\n";
		print WW "newInclusionSeq:\n" . $newInclSeq . "\n";

	}
}
close FF;
close WW;
