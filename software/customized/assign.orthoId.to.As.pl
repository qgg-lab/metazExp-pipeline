#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--taxonId \\\n" .
                "--orthAsConserLevel \\\n" .
                "--orthAsList \\\n" .
		"--checkLineageOrderList \\\n" .
		"--outputAsWithOrthoIdAndConservLevel \n";
	exit;
}

my ($taxonId, $orthAsConserLevel, $orthAsList, $checkLineageOrderList,
$outputAsWithOrthoIdAndConservLevel);

GetOptions(
        'taxonId=s'=>\$taxonId,
	'checkLineageOrderList=s'=>\$checkLineageOrderList,
        'orthAsConserLevel=s'=>\$orthAsConserLevel,
        'orthAsList=s'=>\$orthAsList,
        'outputAsWithOrthoIdAndConservLevel=s'=>\$outputAsWithOrthoIdAndConservLevel,
);

# 将taxonId物种包含的orthoAsId读入hash
# orthAsId 		3702 4577 39947 3847 3880 4081 3711 4113 29760 326968 Angi Mono Poal Dico Capp Sola Faba Rham
# orthA3SS00000012      2    0    0     0    0    0    0    0    0     0      0    0    0    0    0    0    0    0
# orthA3SS00000112      0    0    0     0    0    1    0    1    0     0      0    0    0    0    0    1    0    0
my (%orthoAs, $orthoAsHref, %tmpHash);
$orthoAsHref = \%orthoAs;
my ($line, $i, $tmpName);
my (@nameField, $nameField, @valueField, $valueField);
my ($orthoAsId, $crossSpeciesconserv, $intraspecificConserv, @lineage);
@lineage = split(/,/, $checkLineageOrderList);
open FF, "<$orthAsConserLevel";
$line = <FF>;
chomp($line);
@nameField = split(/\t/, $line);
while($line=<FF>){
	chomp($line);

	@valueField = split(/\t/, $line);
	for($i=0; $i<=$#valueField; $i++){
		$tmpHash{$nameField[$i]} = $valueField[$i];
	}

	# 获得taxonId下的orthoAsId
	next if($tmpHash{$taxonId} == 0);

	$orthoAsId = $tmpHash{"orthAsId"};

	# 扫描conservLevel部分，获得保守水平
	$crossSpeciesconserv = "-";
	for($i=0; $i<=$#lineage; $i++){
		if($tmpHash{$lineage[$i]} == 1){
			$crossSpeciesconserv = $lineage[$i];
			last;
		}
	}
	
	# 获得种内保守性
	$intraspecificConserv = $tmpHash{$taxonId};

	# 登记到hash中
	$orthoAsHref->{$orthoAsId}->{"crossSpeciesconserv"} = $crossSpeciesconserv;
}
close FF;

# 读取orthoAsId和AS之间的映射关系，获得所有AS的orthoAsId，然后在将保守水平添加上去.
# orthAsId        3702    4577    39947   3847    3880    4081    3711    4113    29760   326968
# orthA3SS00000019        -       ZMAYA3SS0000009214,ZMAYA3SS0000001109   -       -       -       -       -       -       -       -
my (@asId, $asId, $asIdList, $innerConservAsIdList);
open WW, ">$outputAsWithOrthoIdAndConservLevel";
print WW join("\t", "ASID", "orthoAsId", "crossSpeciesconserv", "innerConserv") . "\n";
open FF, "<$orthAsList";
$line = <FF>;
chomp($line);
@nameField = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	
	@valueField = split(/\t/, $line);
	for($i=0; $i<=$#valueField; $i++){
		$tmpHash{$nameField[$i]} = $valueField[$i];
	}
	
	# 获得orthoAs名下的AS列表
	$asIdList = $tmpHash{$taxonId};
	next if($asIdList eq "-");

	@asId = ();
	@asId = split(/,/, $asIdList);
	$orthoAsId = $tmpHash{"orthAsId"};

	foreach $asId(@asId){
		# 获得种内保守的AS
		$innerConservAsIdList = "";
		for($i=0; $i<=$#asId; $i++){
			if($asId[$i] ne $asId){
				$innerConservAsIdList .= $asId[$i] . ",";
			}
		}
		if($innerConservAsIdList ne ""){
			$innerConservAsIdList = substr($innerConservAsIdList, 0, length($innerConservAsIdList) - 1);
		}else{
			$innerConservAsIdList = "-";
		}
		if($orthoAsHref->{$orthoAsId}->{"crossSpeciesconserv"} ne "-"){
			print WW join("\t", $asId, $orthoAsId, $orthoAsHref->{$orthoAsId}->{"crossSpeciesconserv"}, $innerConservAsIdList) . "\n";
		}elsif($innerConservAsIdList ne "-"){
			print WW join("\t", $asId, $orthoAsId, "Inner", $innerConservAsIdList) . "\n";
		}
	}
}
close FF;
close WW;
