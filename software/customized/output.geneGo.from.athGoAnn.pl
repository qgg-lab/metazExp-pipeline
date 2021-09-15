#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--geneIdToAthaGeneId \\\n" .
                "--athGoAnno \\\n" .
                "--outputGeneGOAnnoTsv \n";
	exit;
}

my ($geneIdToAthaGeneId, $athGoAnno, $outputGeneGOAnnoTsv);

GetOptions(
        'geneIdToAthaGeneId=s'=>\$geneIdToAthaGeneId,
        'athGoAnno=s'=>\$athGoAnno,
        'outputGeneGOAnnoTsv=s'=>\$outputGeneGOAnnoTsv,
);

my (%geneId2AthGeneId, %athGoList, $line, @field, $geneId, $athGeneId, $goList, $hitGeneId);

open FF, "<$athGoAnno";
while($line=<FF>){
	chomp($line);
	($athGeneId, $goList) = split(/\t/, $line);
	$athGoList{$athGeneId} = $goList;
}
close FF;

open FF, "<$geneIdToAthaGeneId";
open WW, ">$outputGeneGOAnnoTsv";
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	$geneId = $field[0];
	$hitGeneId = $field[1];
	if(exists($athGoList{$hitGeneId})){
		print WW join("\t", $geneId, $athGoList{$hitGeneId}) . "\n";
	}
}
close FF;
close WW;
