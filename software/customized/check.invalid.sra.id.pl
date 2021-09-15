#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--psiDir \\\n" .
                "--inputIdWithoutRltList \\\n" .
		"--outputInvalidIdList \\\n" .
		"--unFinishedIdList \n";
	exit;
}

my ($psiDir, $inputIdWithoutRltList,  $outputInvalidIdList, $unFinishedIdList);

GetOptions(
        'psiDir=s'=>\$psiDir,
        'inputIdWithoutRltList=s'=>\$inputIdWithoutRltList,
        'outputInvalidIdList=s'=>\$outputInvalidIdList,
        'unFinishedIdList=s'=>\$unFinishedIdList,
);

my ($status);
open INVALID, ">$outputInvalidIdList";
open UNFINISHED, ">$unFinishedIdList";
open FF, "<$inputIdWithoutRltList";
while(my $id=<FF>){
	chomp($id);
	$status = "unfinished";
	open TT, "<$psiDir/$id/$id.SeqInfo.txt";
	while(my $line=<TT>){
		if($line=~/\tNA\t\n/){
			$status = "NA";
		}
	}
	close TT;
	if($status eq "unfinished"){
		print UNFINISHED $id . "\n";
	}else{
		print INVALID $id . "\n";
	}
}
close INVALID;
close UNFINISHED;
close FF;
