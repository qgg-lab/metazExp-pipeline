#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--taxonNameAbbrList \\\n" .
                "--dir \n";
	exit;
}

my ($taxonNameAbbrList, $parameterFile, $dir);

GetOptions(
        'taxonNameAbbrList=s'=>\$taxonNameAbbrList,
        'dir=s'=>\$dir,
);

my ($taxonId, $speciesName, $abbr, $line);
open FF, "<$taxonNameAbbrList";
while($line=<FF>){
	chomp($line);
	($speciesName, $abbr, $taxonId) = split(/\t/, $line);
	system("mkdir -p $dir/$taxonId/012-generate-AS-catalog");
	open WW, ">$dir/$taxonId/012-generate-AS-catalog/00000.parameter.of.0000.job1.cfg";
	print WW "export speciesName=\"$speciesName\"\n";
	print WW "export speciesAbbr=$abbr\n";
	close WW;
}
close FF;
