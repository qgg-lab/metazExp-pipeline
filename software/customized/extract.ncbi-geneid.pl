#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--gtfFile \\\n" .
		"--outputNcbiGeneIdTsv \n";
	exit;
}

my ($gtfFile, $outputNcbiGeneIdTsv);

GetOptions(
        'gtfFile=s'=>\$gtfFile,
        'outputNcbiGeneIdTsv=s'=>\$outputNcbiGeneIdTsv,
);

my ($line, @field, @attr, $attr, $geneId, $ncbiGeneId);
my (%geneNcbiGeneId);
open FF, "<$gtfFile";
open WW, ">$outputNcbiGeneIdTsv";
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	next if($field[2] ne "gene");
	
	# db_xref "GeneID:107942799";
	@attr = split(/;/, $field[8]);
	($geneId, $ncbiGeneId) = ("", "");
	foreach $attr(@attr){
		if($attr=~/gene_id "(.*)"/){
			$geneId = $1;
		}
		if($attr=~/db_xref "GeneID:(\d+)"/){
			$ncbiGeneId = $1;
		}
	}
	if($geneId ne "" and $ncbiGeneId ne ""){
		print WW join("\t", $geneId, $ncbiGeneId) . "\n";
	}
}
close FF;
close WW;
