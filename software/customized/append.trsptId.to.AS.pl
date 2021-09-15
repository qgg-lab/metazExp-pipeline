#!/usr/bin/perl
use Getopt::Long;
use strict;
if($#ARGV<0){
	print "\tperl $0 \\\n" . 
		"\t\t--asWithLocalExonCoordsFile total.A5SS.with.exon.coordinates.tsv \\\n" .
		"\t\t--asType A5SS \\\n" .
		"\t\t--gtfFile final.complete.trspt.anno.gtf \\\n" . 
		"\t\t--outputAsWithTrsptIdFile total.A5SS.with.trsptId.tsv\n\n";
	print "This script is used to assign trspts to each AS.\n\n";
	exit;
}

my ($asWithLocalExonCoordsFile, $asType, $gtfFile, $outputAsWithTrsptIdFile);

GetOptions(
	'asWithLocalExonCoordsFile=s'=>\$asWithLocalExonCoordsFile,
	'asType=s'=>\$asType,
	'gtfFile=s'=>\$gtfFile,
	'outputAsWithTrsptIdFile=s'=>\$outputAsWithTrsptIdFile
);

if(uc($asType) ne "A5SS" and uc($asType) ne "A3SS" and uc($asType) ne "RI" and uc($asType) ne "SE" and uc($asType) ne "MXE"){
	print "asType must be specified as A5SS, A3SS, RI, SE or MXE\n";
	exit;
}


my ($line, @field, @tmpArr, %asHash, @fieldName);
my ($chain, $chrId,  $i,  $trspt1, $trspt2, $exon1, $exon2, $exon3, $exon4);
open FF, "<$asFile";
$line=<FF>;
chomp($line);
open WW, ">$outputAsFile";
print WW $line . "\tfirstLocalConcatExons\tsecondLocalConcatExons\n";


