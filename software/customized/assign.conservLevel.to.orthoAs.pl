#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--inputOrthoAsStatisticFile \\\n" .
		"--listLineageTaxonIdList \\\n" .
                "--lineageNameList \\\n" .
		"--outputOrthoAsStatisticFileWithConservLevel \\\n";
	exit;
}

my ($inputOrthoAsStatisticFile, $inputOrthoAsListFile, $lineageNameList, $listLineageTaxonIdList, $outputOrthoAsStatisticFileWithConservLevel, $outputOrthoAsListFileWithConservLevel);

GetOptions(
        'inputOrthoAsStatisticFile=s'=>\$inputOrthoAsStatisticFile,
	'listLineageTaxonIdList=s'=>\$listLineageTaxonIdList,
        'lineageNameList=s'=>\$lineageNameList,
	'outputOrthoAsStatisticFileWithConservLevel=s'=>\$outputOrthoAsStatisticFileWithConservLevel,
);

# 将lineage读入数组
my (@lineageName, $lineageName, @lineageTaxonIdList, $lineageTaxonIdList, @taxonId, $taxonId, $line);
my (%asNum, $asNumHref, @nameField, @valueField, $i, $j,  $outputLine);
$asNumHref = \%asNum;
# Angi Mono Poal Dico Capp Sola Faba 
@lineageName = split(/,/, $lineageNameList);
# 39947_4577_3702_4081_3847 39947_4577 39947_4577 3702_4081_3847 3702 4081 3847 
@lineageTaxonIdList = split(/,/, $listLineageTaxonIdList);

open FF, "<$inputOrthoAsStatisticFile";
$line = <FF>;
chomp($line);
@nameField = split(/\t/, $line);
open WW, ">$outputOrthoAsStatisticFileWithConservLevel";
# orthAsId                3702    3847    39947   4577    4081   Angi Mono Poal Dico Capp Sola Faba
# orthA3SS00000001        0       0       1       0       0      0    0    0    0    0    0    0
print WW join("\t", $line, @lineageName) . "\n";
while($line = <FF>){
	chomp($line);
	@valueField = ();
	@valueField = split(/\t/, $line); 
	
	# 将当前as在各个物种中出现的数量登记到hash中
	%asNum = ();
	for($i=0; $i<=$#valueField; $i++){
		$asNumHref->{$nameField[$i]} = $valueField[$i];
	}

	$outputLine = $line;
	for($i=0; $i<=$#lineageTaxonIdList; $i++){
		@taxonId = ();
		@taxonId = split(/_/, $lineageTaxonIdList[$i]);
		# 检查每个物种taxonId中出现的数量是否都大于0,
		my $checkRlt = 1;
		foreach $taxonId(@taxonId){
			$checkRlt *= $asNumHref->{$taxonId};
		}
		if($checkRlt !=0){
			$outputLine .= "\t1"; 
		}else{
			$outputLine .= "\t0";
		}
	}
	print WW $outputLine . "\n";
}
close FF;
close WW
