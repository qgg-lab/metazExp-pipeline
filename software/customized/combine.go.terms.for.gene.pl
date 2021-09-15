#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--inputGeneGoList \\\n" .
		"--outputCombinedGoList \n";
	exit;
}

my ($inputGeneGoList, $outputCombinedGoList);

GetOptions(
        'inputGeneGoList=s'=>\$inputGeneGoList,
        'outputCombinedGoList=s'=>\$outputCombinedGoList,
);

my (%geneToGoList, $geneId, $go, $line);
open FF, "<$inputGeneGoList";
while($line=<FF>){
	chomp($line);

	($geneId, $go) = split(/\t/, $line);

	if(not exists($geneToGoList{$geneId})){

		$geneToGoList{$geneId} = $go;

	}elsif(index($geneToGoList{$geneId}, $go)<0){

		$geneToGoList{$geneId} .= "," . $go;
	}
}
close FF;

open WW, ">$outputCombinedGoList";
my @geneId = keys(%geneToGoList);
foreach $geneId(@geneId){
	print WW join("\t", $geneId, $geneToGoList{$geneId}) . "\n";
}
close WW;
