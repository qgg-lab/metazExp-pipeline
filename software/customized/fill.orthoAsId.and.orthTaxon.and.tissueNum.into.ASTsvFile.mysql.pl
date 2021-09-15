#!/usr/bin/perl
use strict;
use DBI;
use Getopt::Long;
use DBI;
if($#ARGV<0){
        print "\nperl $0 \\\n" .
		"--asTsvFile origin.as.tsv\\\n" .
		"--orthAsMatrix orthAs.matrix.tsv \\\n" .
		"--psiTsvFile psi.tsv \\\n" .
		"--expInformation experiment.with.tissue.and.prjId.tsv \\\n" .
		"--taxonId 9031\\\n" .
		"--outputAsTsvFile as.with.orthoAs.tsv \n\n";
        exit;
}

my ($asTsvFile, $orthAsMatrix, $taxonId, $expInformation, $psiTsvFile, $outputAsTsvFile);
GetOptions(
	'asTsvFile=s'=>\$asTsvFile,
	'orthAsMatrix=s'=>\$orthAsMatrix,
	'psiTsvFile=s'=>\$psiTsvFile,
	'expInformation=s'=>\$expInformation,
	'taxonId=s'=>\$taxonId,
	'outputAsTsvFile=s'=>\$outputAsTsvFile,
);

my (@fieldTitle, $fieldTitle, @field, $field, $line, @tmp);
my ($orthAsId, $asId, %orthAsTaxon);
my (%asToOrthAs);
open FF, "<$orthAsMatrix";
$line = <FF>;
chomp($line);
@fieldTitle = split(/\t/, $line);
shift(@fieldTitle);
#orthASId        	 9031    		9796    9823    9913    9940
#ORTHA3SSA3SS0000000001  GGALA3SS0000017966     -       -       -       -
while($line = <FF>){
	chomp($line);
	@field = split(/\t/, $line);
	$orthAsId = shift(@field);
	for(my $i=0; $i<=$#field; $i++){
		${$asToOrthAs{$fieldTitle[$i]}}{$field[$i]} = $orthAsId;
		$orthAsTaxon{$orthAsId}.= $fieldTitle[$i] . "," if($field[$i] ne "-");
	}
	$orthAsTaxon{$orthAsId} = substr($orthAsTaxon{$orthAsId}, 0, length($orthAsTaxon{$orthAsId}) -1);
}
close FF;

my ($fieldNameString, $valueString, $expId, $tissue, %expIdToTissue);
# 读取experiment信息，建立expId -> tissue 关系
open FF, "<$expInformation";
while($line = <FF>){
	chomp($line);
	($fieldNameString, $valueString) = split(/_____/, $line);

	@field = split(/, /, $valueString);
	$expId = substr($field[2], 1, length($field[2]) - 2);
	$tissue = substr($field[$#field], 1, length($field[$#field]) - 2);
	$expIdToTissue{$expId} = $tissue;
}
close FF;

# 读取psi，建立asId -> tissueList
my ($inclusion, $skipping, %asToTissueNum);
open FF, "<$psiTsvFile";
while($line = <FF>){
	chomp($line);
	($fieldNameString, $valueString) = split(/_____/, $line);

	@field = split(/, /, $valueString);
	$asId = substr($field[0], 1, length($field[0]) - 2);
	$expId = substr($field[2], 1, length($field[2]) - 2);
	$inclusion = substr($field[4], 1, length($field[4]) - 2);		
	$skipping = substr($field[5], 1, length($field[5]) - 2);
	if($inclusion!=0 or $skipping!=0){
		${$asToTissueNum{$asId}}{$expIdToTissue{$expId}}++;
	}
}
close FF;

#
my ($tissueNumList, @tissue, $tissue, $tissueNum);
open FF, "<$asTsvFile";
open WW, ">$outputAsTsvFile";
while($line=<FF>){
	chomp($line);
	($fieldNameString, $valueString) = split(/_____/, $line);

	@field = split(/, /, $valueString);
	$asId = substr($field[0], 1, length($field[0]) - 2);

	print WW $fieldNameString . ", " . "orthAsId" . ", " . "orthAsTaxon" . ", " . "TissueExpNumList" . ", " . "TissueNum";
	print WW "_____";
	print WW $valueString . ", \"" . ${$asToOrthAs{$taxonId}}{$asId} . "\", \"" . $orthAsTaxon{${$asToOrthAs{$taxonId}}{$asId}} . "\", \"";
	$tissueNumList = "";
	@tissue = ();
	@tissue = keys(%{$asToTissueNum{$asId}});
	$tissueNum = $#tissue + 1;
	foreach $tissue(@tissue){
		$tissueNumList .= $tissue . "(" . ${$asToTissueNum{$asId}}{$tissue} . "),";
	}
	$tissueNumList = substr($tissueNumList, 0, length($tissueNumList)-1);
	print WW $tissueNumList . "\", " . $tissueNum . "\n";
}
close FF;
close WW;
