#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--idList \\\n" .
                "--baseDir \\\n" .
		"--percentCutoff \\\n" .
		"--outputValidIdList \\\n" .
		"--outputInvalidIdList \\\n";
	exit;
}

my ($idList, $baseDir, $percentCutoff, $outputValidIdList, $outputInvalidIdList);

GetOptions(
        'idList=s'=>\$idList,
        'baseDir=s'=>\$baseDir,
	'percentCutoff=s'=>\$percentCutoff,
	'outputValidIdList=s'=>\$outputValidIdList,
	'outputInvalidIdList=s'=>\$outputInvalidIdList
);

my ($alignFile);
my (@id, $id, @line, $line, $totalSpot, $alignPercent);
open FF, "<$idList";
@id=<FF>;
close FF;

open INVALID, ">$outputInvalidIdList";
open VALID, ">$outputValidIdList";
foreach $id(@id){
	chomp($id);
	$alignFile = $baseDir . "/$id/$id.SeqInfo.txt";
	# TotalSpot:29381211      AlignPercent:85.88%
	$totalSpot = 0;
	$alignPercent = 0;
	open FF, "<$alignFile";
	@line=();
	@line=<FF>;
	close FF;
	foreach $line(@line){
		if($line=~/TotalSpot:(\d+)\tAlignPercent:(.*)%\n/){
			$totalSpot = $1;
			$alignPercent = $2;
		}
	}
	if($alignPercent<$percentCutoff){
		print INVALID join("\t", $id) . "\n";
	}else{
		print VALID join("\t", $id) . "\n";
	}
}
close INVALID;
close VALID;
