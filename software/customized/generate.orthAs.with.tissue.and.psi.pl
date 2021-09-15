#!/usr/bin/perl
use strict;
use Getopt::Long;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--tissuePsiFileList 9031.tissue.psi.tsv,9796.tissue.psi.tsv,9823.tissue.psi.tsv,9913.tissue.psi.tsv,9940.tissue.psi.tsv \\\n" .
                "--orthAsMatrix orthAs.matrix.tsv \\\n" .
		"--speciesList \"Gallus gallus,Equus caballus,Sus scrofa,Bos taurus,Ovis aries\"\\\n" .
		"--outputAsTissuePsi orthAs.tissue.psi.tsv \n";
	exit;
}

my ($tissuePsiFileList, $orthAsMatrix, $taxonIdList, $outputAsTissuePsi, $speciesList);
GetOptions(
        'tissuePsiFileList=s'=>\$tissuePsiFileList,
	'orthAsMatrix=s'=>\$orthAsMatrix,
	'speciesList=s'=>\$speciesList,
	'outputAsTissuePsi=s'=>\$outputAsTissuePsi,
       );

my (@tissuePsiFile, $tissuePsiFile, @taxonId, $taxonId, $i, @species, $species);
my (%orthAsToAs, %asToOrthAs, $line);
my (@fieldTitle, $fieldTitle, @field, $field, $orthAsId);
my ($nextFlag, $asId, @asId, @tissue, $tissue);

# 建立as -> orthAs 关系和 orthAs -> as 关系
open FF, "<$orthAsMatrix";
#orthASId        9031    9796    9823    9913    9940
#ORTHA3SSA3SS0000000001  GGALA3SS0000017966      -       -       -       -
$line = <FF>;
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	$orthAsId = shift(@field);

	# 检查是否需要登记该orthAs
	$nextFlag = 0;
	foreach $field(@field){
		$nextFlag++ if($field ne "-");
	}
	next if($nextFlag <= 1);

	# 存在两个以上不同物种的as在该该orthAs中时，给予登记
	for($i=0; $i<=$#field; $i++){
		if($field[$i] ne "-"){
			$asToOrthAs{$field[$i]} = $orthAsId;
		}
	}
}
close FF;

# 将不同组织、不同物种的psi组织到psi哈希表中
my (%psi, @tissuePsi, $tissuePsi, $tissue, @psi, $psi); 
@tissuePsiFile = split(/,/, $tissuePsiFileList);
@species = split(/,/, $speciesList);

for($i=0; $i<=$#species; $i++){
	$tissuePsiFile = $tissuePsiFile[$i];
	$species = $species[$i];
	open FF, "<$tissuePsiFile";
	while($line = <FF>){
		#GGALA3SS0000009692  blastoderm,0.38205499276411,0.0717391304347826 liver,0.983,..
		chomp($line);
		@tissuePsi = ();
		@tissuePsi = split(/\t/, $line);
		$asId = shift(@tissuePsi);
		next if(not exists($asToOrthAs{$asId}));
		$orthAsId = $asToOrthAs{$asId};
		foreach $tissuePsi(@tissuePsi){
			@psi = ();
			@psi = split(/,/, $tissuePsi);
			$tissue = shift(@psi);
			${$psi{$orthAsId}}{$tissue} .= $species . "," . join(",", @psi) . "#";
		}
	}
	close FF;
}

open WW, ">$outputAsTissuePsi";
my (@tissueSpeciesPsi, $tissueSpeciesPsi);
my (@speciesPsi, $speciesPsi);
my @orthAsId = keys(%psi);
foreach $orthAsId(@orthAsId){

	$nextFlag = 0;
	
	# 检查该orthAs是否存在这样的tissue：有至少2个物种。如果不存在那么放弃该orthAs。
	# 如果存在，那么只输出该orthAs的包含2个以上物种tissue的psi
	# 输出的格式如下：
	# 组织间用\t分隔；组织内物种间用#分隔
	# orthAsId \t liver:cattle,0.9,0.3,0.1#sheep,0.1,0.4,0.7#chicken,0.1,0.3,0.7 \t lung:cattle,0.1,0.2,0.3#sheep,0.3,0.5,0.2
	@tissue = ();
	@tissue = keys(%{$psi{$orthAsId}});
	foreach $tissue(@tissue){
		@species = ();
		@species = split(/#/, ${$psi{$orthAsId}}{$tissue});
		$nextFlag = 1 if($#species>=1);
	}

	next if($nextFlag == 0);

	print WW $orthAsId;
	@tissue = ();
	@tissue = keys(%{$psi{$orthAsId}});
	foreach $tissue(@tissue){
		@species = ();
		@species = split(/#/, ${$psi{$orthAsId}}{$tissue});

		next if($#species<=0);

		print WW "\t" . $tissue . ":";
		print WW ${$psi{$orthAsId}}{$tissue};	
	}
	print WW "\n";
}
close WW;
