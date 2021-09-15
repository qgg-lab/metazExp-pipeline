#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
        print "\nperl $0 \\\n" .
                "--sqlFile 9925.speciesTable.sql \\\n" .
                "--taxonId 9925 \\\n" . 
		"--tmpDir ./tmpDir \n\n\n";
	exit;
}
my ($sqlFile, $taxonId, $tmpDir, $line, @valueItems, $totalValueItems, $valueItem, $prefixString);
my (@fields);
GetOptions(
        'sqlFile=s'=>\$sqlFile,
        'taxonId=s'=>\$taxonId,
	'tmpDir=s'=>\$tmpDir,
);

open FF, "<$sqlFile";
system("mkdir -p " . $tmpDir);
open WW, ">$tmpDir/tmp.$taxonId.speciesTable.sql";
while($line=<FF>){
	if($line=~/^DROP TABLE IF EXISTS.*;/){
		print WW "/*" . $line;
	}elsif($line=~/^CREATE TABLE/){
		print WW "/*" . $line;
	}elsif($line=~/^(INSERT INTO .*? VALUES \()(.*)\);\n/){
		$prefixString = $1;
		$totalValueItems = $2;
		@valueItems = split(/\),\(/, $totalValueItems);
		foreach $valueItem(@valueItems){
			@fields = ();
			@fields = split(/,/, $valueItem);
			if($fields[1]=~/'$taxonId'/){
				print WW $prefixString . $valueItem . ");\n";
			}
		}
	}else{
		print WW $line;
	}
}
close FF;
system("rm -rf $sqlFile");
system("mv " . "$tmpDir/tmp.$taxonId.speciesTable.sql " . $sqlFile);
