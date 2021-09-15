#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--geneSpecificExprTpmTsv \\\n" .
                "--outputGeneSpecificExprTsv \n";
	exit;
}

my ($geneSpecificExprTpmTsv, $outputGeneSpecificExprTsv);

GetOptions(
        'geneSpecificExprTpmTsv=s'=>\$geneSpecificExprTpmTsv,
        'outputGeneSpecificExprTsv=s'=>\$outputGeneSpecificExprTsv,
);

my ($avgValueList);
my (@field, $line, $GeneId);
my ($SpecificTissue, $SpecificType, $SpecificPvalue, $SpecificValueList, $SpecificAvgValue, $OtherTissueValueInfoList);
my ($fieldNameList, $valueList, $tissueNum);

open WW, ">$outputGeneSpecificExprTsv";
open FF, "<$geneSpecificExprTpmTsv";
#pvalue0.01Rlt pvalue0.03Rlt pvalue0.05Rlt geneId High/Low tissue avgExpTpm|TpmList 1stTissue|pvalue|avg(TPM)|TPMList 2ndTissue|pvalue|avg(TPM)|TPMList
<FF>;
while($line=<FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);

	next if($field[0] eq "N" and $field[1] eq "N" and $field[2] eq "N");

	$SpecificPvalue = "0.05" if($field[2] eq "Y");
	$SpecificPvalue = "0.03" if($field[1] eq "Y");
	$SpecificPvalue = "0.01" if($field[0] eq "Y");

	shift(@field);
	shift(@field);
	shift(@field);
	$GeneId	= shift(@field);
	$SpecificType = shift(@field);	
	$SpecificTissue = shift(@field);
	$avgValueList = shift(@field);
	($SpecificAvgValue, $SpecificValueList) = split(/\|/, $avgValueList);
	$tissueNum = $#field+1;
	$OtherTissueValueInfoList = join(";", @field);

	$fieldNameList = join("###", "TrsptId", "SpecificTissue", "SpecificType", "SpecificPvalue", "SpecificAvgValue", "SpecificValueList", "OtherTissueValueInfoList", "tissueNum");
	$valueList = join("###", $GeneId, $SpecificTissue, $SpecificType, $SpecificPvalue, $SpecificAvgValue, $SpecificValueList, $OtherTissueValueInfoList, $tissueNum);
	print WW $fieldNameList . "___" . $valueList . "\n";
}
close FF;
close WW;
