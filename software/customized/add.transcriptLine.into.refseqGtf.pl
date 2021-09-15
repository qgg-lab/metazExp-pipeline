#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--refseqGtf \\\n" .
                "--outputGtf \n";

	exit;
}

my ($refseqGtf, $outputGtf);

GetOptions(
        'refseqGtf=s'=>\$refseqGtf,
        'outputGtf=s'=>\$outputGtf,
);

my (%trsptSpan);
