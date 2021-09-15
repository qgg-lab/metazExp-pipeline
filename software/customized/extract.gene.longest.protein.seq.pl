#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--longestProteome \\\n" .
                "--geneIdList \\\n" .
                "--outputSeqFile \n";
	exit;
}

my ($longestProteome, $geneIdList, $outputSeqFile);

GetOptions(
        'longestProteome=s'=>\$longestProteome,
        'geneIdList=s'=>\$geneIdList,
        'outputSeqFile=s'=>\$outputSeqFile,
);


