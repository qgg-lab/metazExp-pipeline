#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--orthoGeneGroupTsv \\\n" .
                "--taxonIdList \\\n" .
		"--speciesList \\\n" .
                "--speciesAbbrList \\\n" .
		"--orderList \\\n" .
		"--outputMysqlTsv \n";
	exit;
}

my ($orthoGeneGroupTsv, $taxonIdList, $speciesList, $speciesAbbrList, $orderList, $outputMysqlTsv);

GetOptions(
        'orthoGeneGroupTsv=s'=>\$orthoGeneGroupTsv,
        'taxonIdList=s'=>\$taxonIdList,
	'speciesList=s'=>\$speciesList,
        'speciesAbbrList=s'=>\$speciesAbbrList,
        'orderList=s'=>\$orderList,
	'outputMysqlTsv=s'=>\$outputMysqlTsv,
);
my (@taxonId, $taxonId, @speciesAbbr, $speciesAbbr, @order, $order, $line, $orthoGeneGroupId, $modifiedGeneIdList, @modifiedGeneId, $modifiedGeneId, $geneId, $i, @species, $species);

@taxonId = split(/,/, $taxonIdList);
@speciesAbbr = split(/,/, $speciesAbbrList);
@species = split(/,/, $speciesList);
@order = split(/,/, $orderList);

# 建立：ABBR->TAXONID;ABBR->ORDER的映射hash
my (%abbrToTaxonIdToOrder, $abbrToTaxonIdToOrderHref);
$abbrToTaxonIdToOrderHref = \%abbrToTaxonIdToOrder;
for($i=0; $i<=$#speciesAbbr; $i++){
	$abbrToTaxonIdToOrderHref->{$speciesAbbr[$i]}->{"taxonId"}=$taxonId[$i];
	$abbrToTaxonIdToOrderHref->{$speciesAbbr[$i]}->{"order"}=$order[$i];
	$abbrToTaxonIdToOrderHref->{$speciesAbbr[$i]}->{"species"}=$species[$i];
}


my (@nameField, @valueField,  $valueString, %hashTmp, $tmpHref, $abbr, $taxonId);
$tmpHref=\%hashTmp;

my $nameString=join(", ", "geneId", "orthoGeneGroupId", "taxonId", "species", "speciesAbbr", "orderName");

open WW, ">$outputMysqlTsv";


# 打开直系同源groupTsv
open FF, "<$orthoGeneGroupTsv";

# 获得每一列的名称
# Orthogroup   AHAL    ZJUJ   AHYP    ALYR    ATHA    BNAP 
$line=<FF>;
chomp($line);
@nameField = split(/\t/, $line);

while($line=<FF>){
	
	# OG0039044 ZJUJ_LOC112492972_pep001, ZJUJ_LOC112492974_pep001    ATHA_xxx_pep001
	chomp($line);

	@valueField = split(/\t/, $line);

	%hashTmp = ();

	# 把行中每一列值放到hash中
	for($i=0; $i<=$#nameField; $i++){
		my $abbr = $nameField[$i];
		my $geneIdList= $valueField[$i];
		$tmpHref->{$abbr} = $geneIdList;
		print "abbr:" . $abbr . ":\n";
		print "geneIdList:" . $geneIdList . "\n";
		print "comOut:$abbr:" . $tmpHref->{$abbr} . "\n";
		<STDIN>;
	}
	
	$orthoGeneGroupId = $tmpHref->{"Orthogroup"};
	
	# 扫描每个物种
	for($i=0; $i<=$#speciesAbbr; $i++){
		# 获得该物种的缩写、分类号和目名称
		$abbr = $speciesAbbr[$i];

		$taxonId = $abbrToTaxonIdToOrderHref->{$abbr}->{"taxonId"};
		$order = $abbrToTaxonIdToOrderHref->{$abbr}->{"order"};
		$species = $abbrToTaxonIdToOrderHref->{$abbr}->{"species"};

#		print join("\t", $abbr, $taxonId, $order, $tmpHash{$abbr});
#		<STDIN>;
		next if($tmpHref->{$abbr} eq "");

		# 获得该物种下所有基因的ID
		$modifiedGeneIdList = $tmpHref->{$abbr};
		@modifiedGeneId = ();
		@modifiedGeneId = split(/, /, $modifiedGeneIdList);

		foreach $modifiedGeneId(@modifiedGeneId){
			#print $modifiedGeneId;
			#<STDIN>;
			if($modifiedGeneId=~/$abbr(_)(.*?)_pep\d+/){
				$geneId = $2;
				$valueString = join(", ", $geneId, $orthoGeneGroupId, $taxonId, $species, $abbr, $order);
				print WW $nameString . "___" . $valueString . "\n";
			}						
		}
	}
}
close FF;
close WW;
