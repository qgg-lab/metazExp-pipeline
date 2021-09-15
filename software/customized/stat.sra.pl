#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--filterExptInfoTsv \\\n" .
                "--outputStatTsv \n";
	exit;
}

my ($filterExptInfoTsv, $outputStatTsv);

GetOptions(
        'filterExptInfoTsv=s'=>\$filterExptInfoTsv,
        'outputStatTsv=s'=>\$outputStatTsv,
);

# Ecotype Cultivar        Genotype        Tissue  SubTissue       Development     Treatment       treatmentGroup  Experiment      Study   dataSource      Base    Layout  SpotsNum        ReadNum SpotLen ReadLen Gather  AS      Assemble        RunList Phenotype       alignPercent    mappedBases     mappedReadNum   detectedReadLen libraryType     phredScore
#
my (%expt, %tissue, %study, @run, $runNum, $totalBase, $mappedBase, $line, @field, $field, @value, $value, %tmpExpt);
open FF, "<$filterExptInfoTsv";
$line=<FF>;
chomp($line);
@field = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@value = split(/\t/, $line);
	for(my $i=0; $i<=$#value; $i++){
		$tmpExpt{$field[$i]} = $value[$i];
	}
	$study{$tmpExpt{"Study"}} = 1;
	$expt{$tmpExpt{"Experiment"}} = 1;
	# run num
	@run = ();
	@run = split(/,/, $tmpExpt{"RunList"});
	$runNum += $#run + 1;
	
	# tissue
	$tissue{$tmpExpt{"Tissue"}} = 1;

	# totalBase
	$totalBase += $tmpExpt{"Base"};

	# mappedBase
	$mappedBase += $tmpExpt{"mappedBases"};
	
	
}
close FF;

open WW, ">$outputStatTsv";
my (@study, @expt, @tissue);
@study = keys(%study);
@expt = keys(%expt);
@tissue = keys(%tissue);
print WW join("\t", "Study", "Experiment", "Run", "Total Bases", "Mapped Bases") . "\n";
print WW join("\t", $#study+1, $#expt+1, $runNum, sprintf("%.1f",$totalBase), sprintf("%.1f",$mappedBase)) . "\n";
close WW;
