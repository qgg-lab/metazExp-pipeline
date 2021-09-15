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
	'outputValidIdList=s'=>\$outputValidIdList,
	'outputInvalidIdList=s'=>\$outputInvalidIdList
);

my ($alignFile);
my (@id, $id, @line, $line, $naFlag);
open FF, "<$idList";
@id=<FF>;
close FF;

open INVALID, ">$outputInvalidIdList";
open VALID, ">$outputValidIdList";
foreach $id(@id){
	chomp($id);
	$alignFile = $baseDir . "/$id/$id.SeqInfo.txt";
	# ExpId   RunId   Layout  ReadLen PhredScore      PhredDetail     Library LibraryDetail
	# SRX1485840      SRR3019331      PAIRED  75      33      NegativeScoreBaseCount(%):47729(6.43)   NA 
	$naFlag = 0;
	open FF, "<$alignFile";
	@line=();
	@line=<FF>;
	close FF;
	foreach $line(@line){
		if($line=~/\tNA\t/){
			$naFlag = 1;
		}
	}
	if($naFlag==1){
		print INVALID join("\t", $id) . "\n";
	}else{
		print VALID join("\t", $id) . "\n";
	}
}
close INVALID;
close VALID;
