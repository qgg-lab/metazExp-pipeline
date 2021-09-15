#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--specialAsTsvList tissueSpecialHighPsi.tsv,tissueSpeciesLowPsi.tsv \\\n" .
		"--directionList High,Low\\\n" .
		"--experimentAttr  Tissue\\\n" .
		"--outputAsMysqlTsv 3702.tissue.special.as.mysql.insert.tsv\n";
	exit;
}

my ($specialAsTsvList, $directionList, $experimentAttr, $outputAsMysqlTsv);

GetOptions(
	'specialAsTsvList=s'=>\$specialAsTsvList,
	'directionList=s'=>\$directionList,
	'experimentAttr=s'=>\$experimentAttr,
	'outputAsMysqlTsv=s'=>\$outputAsMysqlTsv,
);

my ($line);

open WW, ">$outputAsMysqlTsv";
# 读取组织特异的AS，将组织特异性赋给AS
# ASID    Tissue=MaxPsi   AllFieldTypePsiList
# ATHARI0000005333        meristem=0.656  seedling=0.027  stem=0.022 ...
my ($nameString, $valueString, $i, $asId);
my (@specialAsTsv, $specialAsTsv, @typePsi, $typePsi, $type, $psi, @direction, $direction);

@specialAsTsv = split(/,/, $specialAsTsvList);
@direction=split(/,/, $directionList);
for($i=0; $i<=$#direction; $i++){
	$specialAsTsv = $specialAsTsv[$i];
	$direction = $direction[$i];

	$nameString = join(", ", "asId", $experimentAttr, "specialTag", "psiValue");
	open FF, "<$specialAsTsv";
	<FF>;
	while($line=<FF>){
		chomp($line);
		@typePsi = ();
		@typePsi = split(/\t/, $line);
		$asId = $typePsi[0];
		# 第1个为特异的psi
		($type, $psi) = split(/=/, $typePsi[1]);
		$valueString = join(", ", $asId, $type, $direction, $psi);
		$type=~tr/ /_/;
		print WW join("___", $nameString, $valueString) . "\n";
		for(my $j=2; $j<=$#typePsi; $j++){
			($type, $psi) = split(/=/, $typePsi[$j]);
			$type=~tr/ /_/;
			$valueString = join(", ", $asId, $type, "detected", $psi);
			print WW join("___", $nameString, $valueString) . "\n";
		}

	}
	close FF;
}

close WW;
