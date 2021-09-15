#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--gtfFile \\\n" .
		"--outputProtIdList \\\n" .
		"--outputEntrezIdTsv \n";
	exit;
}

my ($gtfFile, $outputEntrezIdTsv, $outputProtIdList);

GetOptions(
        'gtfFile=s'=>\$gtfFile,
	'outputProtIdList=s'=>\$outputProtIdList,
        'outputEntrezIdTsv=s'=>\$outputEntrezIdTsv,
);

my ($line, @field, @attr, $attr, $geneId, $entrezId, $protId);
my (%geneId2EntrezId, %protId2GeneId);
open FF, "<$gtfFile";
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	
	# db_xref "GeneID:107942799";
	@attr = split(/;/, $field[8]);
	($geneId, $entrezId) = ("", "");
	foreach $attr(@attr){
		if($attr=~/gene_id "(.*)"/){
			$geneId = $1;
		}
		if($attr=~/db_xref "GeneID:(\d+)"/){
			$entrezId = $1;
		}
		if($attr=~/protein_id "(.*)"/){
			$protId = $1;
		}
	}
	if($geneId ne "" and $entrezId ne ""){
		$geneId2EntrezId{$geneId} = $entrezId;
	}

	if($geneId ne "" and $protId ne ""){
		$protId2GeneId{$protId} = $geneId;
	}
}
close FF;


# 输出ncbi中的geneId到entrezId之间的映射
open EE, ">$outputEntrezIdTsv";
my @geneId = keys(%geneId2EntrezId);
foreach my $geneId(@geneId){
	print EE join("\t", $geneId, $geneId2EntrezId{$geneId}) . "\n";
}
close EE;

open WW, ">$outputProtIdList";
my @protId = keys(%protId2GeneId);
foreach $protId(@protId){
	print WW join("\t", $protId, $protId2GeneId{$protId}) . "\n";
}
close WW;

