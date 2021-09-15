#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--asSpecificPsiTsv \\\n" .
                "--outputAsSpecificPsiTsv \n";
	exit;
}

my ($asSpecificPsiTsv, $outputAsSpecificPsiTsv);

GetOptions(
        'asSpecificPsiTsv=s'=>\$asSpecificPsiTsv,
        'outputAsSpecificPsiTsv=s'=>\$outputAsSpecificPsiTsv,
);

my ($avgValueList);
my (@field, $line, $asId);
my ($SpecificTissue, $SpecificType, $SpecificPvalue, $SpecificValueList, $SpecificAvgValue, $OtherTissueValueInfoList, $specificTissueInfo);
my ($fieldNameList, $valueList, $tissueNum);

open WW, ">$outputAsSpecificPsiTsv";
open FF, "<$asSpecificPsiTsv";
# pvalue0.01 pvalue0.03 pvalue0.03 ASID High/Low Tissue=AVG(PSI)|PSIList 1stTissue|pvalue|avg(PSI)|PSIlist
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
	$asId	= shift(@field);
	$SpecificType = shift(@field);	

	$specificTissueInfo = shift(@field);
	$specificTissueInfo=~/(.*)=(.*)\|(.*)/;
	($SpecificTissue, $SpecificAvgValue, $SpecificValueList) = ($1, $2, $3);

	$tissueNum = $#field + 1;
	$OtherTissueValueInfoList = join(";", @field);


	$fieldNameList = join("###", "ASID", "SpecificTissue", "SpecificType", "SpecificPvalue", "SpecificAvgValue", "SpecificValueList", "OtherTissueValueInfoList", "tissueNum");
	$valueList = join("###", $asId, $SpecificTissue, $SpecificType, $SpecificPvalue, $SpecificAvgValue, $SpecificValueList, $OtherTissueValueInfoList, $tissueNum);
	print WW $fieldNameList . "___" . $valueList . "\n";
}
close FF;
close WW;

