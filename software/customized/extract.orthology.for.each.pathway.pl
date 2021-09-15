#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--keggpathwayMysqlTsv \\\n" .
                "--outputOrthologList \n";
	exit;
}

my ($keggpathwayMysqlTsv, $outputOrthologList);

GetOptions(
        'keggpathwayMysqlTsv=s'=>\$keggpathwayMysqlTsv,
        'outputOrthologList=s'=>\$outputOrthologList,
);

my (%pathwayIdToOrthologyList, $nameList, @name, $name, $valueList, @value, $value, %hash);

open FF, "<$keggpathwayMysqlTsv";
while(my $line=<FF>){
	chomp($line);
	($nameList, $valueList) = split(/___/, $line);
	@name = split(/\t/, $nameList);
	@value = split(/\t/, $valueList);
	for(my $i=0; $i<=$#name; $i++){
		$hash{$name[$i]} = $value[$i];
	}
	if(index($pathwayIdToOrthologyList{$hash{"pathwayId"}}, $hash{"orthology"})>=0){
		#$pathwayIdToOrthologyList{$hash{"pathwayId"}} .= $hash{"orthology"} . "\t";
	}else{
		$pathwayIdToOrthologyList{$hash{"pathwayId"}} .= $hash{"orthology"} . "\t";
	}
}
close FF;

open WW, ">$outputOrthologList";
my @pathwayId = keys(%pathwayIdToOrthologyList);
foreach my $pathwayId(@pathwayId){
	print WW join("\t", $pathwayId, $pathwayIdToOrthologyList{$pathwayId}) . "\n";
}
close WW;
