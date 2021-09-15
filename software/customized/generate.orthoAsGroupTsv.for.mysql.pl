#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--orthoAsGroupTsv \\\n" .
                "--taxonIdList \\\n" .
		"--speciesList \\\n" .
                "--speciesAbbrList \\\n" .
		"--orderList \\\n" .
		"--outputMysqlTsv \\\n" .
		"--outputSpeciesNumInOrthoGroup \n";
	exit;
}

my ($orthoAsGroupTsv, $taxonIdList, $speciesList, $speciesAbbrList, $orderList, $outputMysqlTsv, $outputSpeciesNumInOrthoGroup);

GetOptions(
        'orthoAsGroupTsv=s'=>\$orthoAsGroupTsv,
        'taxonIdList=s'=>\$taxonIdList,
	'speciesList=s'=>\$speciesList,
        'speciesAbbrList=s'=>\$speciesAbbrList,
        'orderList=s'=>\$orderList,
	'outputMysqlTsv=s'=>\$outputMysqlTsv,
	'outputSpeciesNumInOrthoGroup=s'=>\$outputSpeciesNumInOrthoGroup,
);
my (@taxonId, $taxonId, @speciesAbbr, $speciesAbbr, @order, $order, $line, $orthoAsGroupId, $asIdList, @asId, $asId, $i, @species, $species);

@taxonId = split(/,/, $taxonIdList);
@speciesAbbr = split(/,/, $speciesAbbrList);
@species = split(/,/, $speciesList);
@order = split(/,/, $orderList);

# 建立：ABBR->TAXONID;ABBR->ORDER的映射hash
my (%taxonIdToAbbrToOrder, $taxonIdToAbbrToOrderHref);
$taxonIdToAbbrToOrderHref = \%taxonIdToAbbrToOrder;
for($i=0; $i<=$#taxonId; $i++){
	$taxonIdToAbbrToOrderHref->{$taxonId[$i]}->{"speciesAbbr"}=$speciesAbbr[$i];
	$taxonIdToAbbrToOrderHref->{$taxonId[$i]}->{"order"}=$order[$i];
	$taxonIdToAbbrToOrderHref->{$taxonId[$i]}->{"species"}=$species[$i];
}


# xxx
my (@nameField, @valueField,  $valueString, %tmpHash, $taxonId, $speciesNum);
my $nameString=join(", ", "asId", "orthoAsGroupId", "taxonId", "species", "speciesAbbr", "orderName", "speciesNum");
open WW, ">$outputMysqlTsv";
open NUM, ">$outputSpeciesNumInOrthoGroup";
# 打开直系同源groupTsv
open FF, "<$orthoAsGroupTsv";
# orthAsId        108875  13443   225117  2711    28526   29729   29730   29760    
$line=<FF>;
chomp($line);
@nameField = split(/\t/, $line);
while($line=<FF>){
# orthA3SS00000049        -    NTABA3SS0000009859,NTABA3SS0000012413   STUBA3SS0000003168 
	chomp($line);
	@valueField = split(/\t/, $line);

	# 统计每个orthoGroup中物种数量
	$speciesNum = &getSpeciesNum($line);

	for($i=0; $i<=$#valueField; $i++){
		$tmpHash{$nameField[$i]} = $valueField[$i];
	}

	$orthoAsGroupId = $tmpHash{"orthAsId"};
	
	# 扫描每个物种
	for($i=0; $i<=$#taxonId; $i++){
		# 获得该物种的缩写、分类号和目名称
		$taxonId = $taxonId[$i];
		$speciesAbbr = $taxonIdToAbbrToOrderHref->{$taxonId}->{"speciesAbbr"};
		$order = $taxonIdToAbbrToOrderHref->{$taxonId}->{"order"};
		$species = $taxonIdToAbbrToOrderHref->{$taxonId}->{"species"};
		next if($tmpHash{$taxonId} eq "-");

		# 获得该物种下所有AS的ID
		$asIdList = $tmpHash{$taxonId[$i]};
		@asId = ();
		@asId = split(/,/, $asIdList);
		foreach $asId(@asId){
			$valueString = join(", ", $asId, $orthoAsGroupId, $taxonId, $species, $speciesAbbr, $order, $speciesNum);
			print WW $nameString . "___" . $valueString . "\n";
		}
	}
}
close FF;
close WW;
close NUM;

sub getSpeciesNum{
	my ($line) = @_;
	my (@tt, $i, $speciesNum);
	$speciesNum = 0;
	@tt = split(/\t/, $line);
	for($i=1; $i<=$#tt; $i++){
		if($tt[$i] ne "-"){
			$speciesNum++;
		}
	}
	return $speciesNum;
}
