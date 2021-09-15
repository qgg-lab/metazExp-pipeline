#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--duplicationGroup \\\n" .
		"--speciesAbbr Atha\\\n" .
                "--duplicationGroupPrefix DPG\\\n" .
                "--outputDG \n";
	exit;
}

my ($duplicationGroup, $duplicationGroupPrefix, $speciesAbbr, $outputDG);

GetOptions(
        'duplicationGroup=s'=>\$duplicationGroup,
        'duplicationGroupPrefix=s'=>\$duplicationGroupPrefix,
	'speciesAbbr=s'=>\$speciesAbbr,
        'outputDG=s'=>\$outputDG,
);

my ($line, @complexGeneId, $complexGeneId, %origGeneId, @origGeneId, $origGeneId, $complexGeneNum, $origGeneNum, $dupGroupNum, $origGeneIdList);
open WW, ">$outputDG";
print WW join("\t", "duplicationGroupId", "geneList", "complextGeneNum", "origGeneNum") . "\n";
open FF, "<$duplicationGroup";
#Atha_AT3G61790_pep002, Atha_AT3G61800.up.cluster1_pep001, Atha_AT4G27880_pep002
while($line=<FF>){
	chomp($line);
	%origGeneId = ();
	@complexGeneId = ();
	@complexGeneId = split(/, /, $line);
	$complexGeneNum = $#complexGeneId + 1;
	foreach $complexGeneId(@complexGeneId){
		if($complexGeneId=~/($speciesAbbr)_(.*)_pep\d+/){
			$origGeneId = $2;
			$origGeneId{$origGeneId} = 1;
		}
	}

	# 获得原始基因的编号
	@origGeneId = keys(%origGeneId);
	$origGeneNum = $#origGeneId + 1;
	
	# 输出原始基因编号
	$origGeneIdList = join(",", @origGeneId);
	$dupGroupNum++;
	print WW join("\t", $duplicationGroupPrefix . sprintf("%05d", $dupGroupNum), $origGeneIdList, $complexGeneNum, $origGeneNum) . "\n";
}
close FF;
