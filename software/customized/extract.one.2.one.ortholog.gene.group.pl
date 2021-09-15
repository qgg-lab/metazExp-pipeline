#!/usr/bin/perl
my $orthologTb = $ARGV[0];
my $abbrList = $ARGV[1];
my $outputOneByOne = $ARGV[2];
my (@abbr, $abbr);

if($#ARGV<2){
	print $0 . " Orthogroups.GeneCount.tsv ATHA,GMAX,NATT one.2.one.group.list\n";
}
@abbr = split(/,/, $abbrList);

open WW, ">$outputOneByOne";
my (@nameField, @valueField, $line, %tmpHash,  $checkRlt);
open FF, "<$orthologTb";
$line = <FF>;
chomp($line);
@nameField = split(/\t/, $line);
# Orthogroup      AHAL    AHYP    ALYR    ATHA
# OG0000000       126     424     19      14
while($line=<FF>){
	chomp($line);
	@valueField = split(/\t/, $line);
	for($i=0; $i<=$#valueField; $i++){
		$tmpHash{$nameField[$i]} = $valueField[$i];
	}

	# 检查关注的物种中基因数量是否都为1
	$checkRlt = 1;
	foreach $abbr(@abbr){
		if($tmpHash{$abbr}!=1){
			$checkRlt = 0;
		}
	}

	if($checkRlt == 1){
		print WW $tmpHash{"Orthogroup"} . "\n";
	}
}
close FF;
close WW;
