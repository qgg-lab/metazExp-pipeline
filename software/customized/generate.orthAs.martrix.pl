#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--orthAsMatrixFileList vertebrate.orth.as.tsv,mammal.orth.as.tsv,artiodactyla.orth.as.tsv,bovidae.orth.as.tsv \\\n" .
		"--lineageLabelList VER,MAM,ART,BOV\\\n" .
		"--outputOrthAsMatrix orthAs.matrix.tsv \n";
	exit;
}

my ($orthAsMatrixFileList, $lineageLabelList, $outputOrthAsMatrix);
GetOptions(
        'orthAsMatrixFileList=s'=>\$orthAsMatrixFileList,
        'lineageLabelList=s'=>\$lineageLabelList,
        'outputOrthAsMatrix=s'=>\$outputOrthAsMatrix,
);


my (@orthAsMatrixFile, $orthAsMatrixFile, @lineageLabel, $lineageLabel, %taxon, @taxon, $taxon);
@orthAsMatrixFile = split(/,/, $orthAsMatrixFileList);
@lineageLabel = split(/,/, $lineageLabelList);

my ($orthAsMatrixFile, %orthAs);
my (@fieldTitle, $fieldTitle, @field, $field, $line, $originOrthAsId, $newOrthAsId, $checkFlag, %orthAsRegister);
my ($orthAsNum);
my ($i, $j);

for($j=0; $j<=$#lineageLabel; $j++){
	$orthAsMatrixFile = $orthAsMatrixFile[$j];
	$lineageLabel = $lineageLabel[$j];
	
	open FF, "<$orthAsMatrixFile";
	# orthASId        9031    9796    9823    9913    9940
	# ORTHA3SSA3SS0000000001  GGALA3SS0000017966      -       -       -       - 
	$line = <FF>;
	chomp($line);
	@fieldTitle = ();
	@fieldTitle = split(/\t/, $line);
	shift(@fieldTitle);
	
	foreach $field(@fieldTitle){
		$taxon{$field} = 1 if(not(exists($taxon{$field})));
	}

	while($line=<FF>){
		chomp($line);
		@field = ();
		@field = split(/\t/, $line);
		$originOrthAsId = shift(@field);

		$checkFlag = 1;
		foreach $field(@field){
			$checkFlag = 0 if($field eq "-");
		}
		next if($checkFlag == 0);

		$orthAsNum++;
		if($originOrthAsId=~/ORTHA3SS/){
			$newOrthAsId = $lineageLabel . "A3SS" . sprintf("%08d", $orthAsNum);
		}elsif($originOrthAsId=~/ORTHA5SS/){
			$newOrthAsId = $lineageLabel . "A5SS" . sprintf("%08d", $orthAsNum);
		}elsif($originOrthAsId=~/ORTHSE/){
			$newOrthAsId = $lineageLabel . "SE" . sprintf("%08d", $orthAsNum);
		}elsif($originOrthAsId=~/ORTHRI/){
			$newOrthAsId = $lineageLabel . "RI" . sprintf("%08d", $orthAsNum);
		}elsif($originOrthAsId=~/ORTHMXE/){
			$newOrthAsId = $lineageLabel . "MXE" . sprintf("%08d", $orthAsNum);
		}
	
		for(my $i=0; $i<=$#fieldTitle; $i++){
			if(not(exists($orthAsRegister{$field[$i]}))){
				${$orthAs{$newOrthAsId}}{$fieldTitle[$i]} = $field[$i];
				$orthAsRegister{$field[$i]} = 1;
			}
		}
	}
	close FF;
}


open WW, ">$outputOrthAsMatrix";
my $orthAsId;
my @orthAsId = keys(%orthAs);
@taxon = keys(%taxon);
@taxon = sort { $a cmp $b } @taxon;

print WW "orthAsId";
foreach $taxon(@taxon){
	print WW "\t" .  $taxon;
}
print WW "\n";


foreach $orthAsId(@orthAsId){
	print WW $orthAsId;
	foreach $taxon(@taxon){
		if(exists(${$orthAs{$orthAsId}}{$taxon})){
			print WW "\t" . ${$orthAs{$orthAsId}}{$taxon};
		}else{
			print WW "\t-";
		}		
	}
	print WW "\n";
}
close WW;
