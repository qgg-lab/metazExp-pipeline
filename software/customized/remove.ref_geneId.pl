#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--origTranscriptome \\\n" .
		"--newTranscriptome \n";
	exit;
}

my ($origTranscriptome, $newTranscriptome);
my ($line, @field, @attr, $attrString, $attr);
GetOptions(
        'origTranscriptome=s'=>\$origTranscriptome,
        'newTranscriptome=s'=>\$newTranscriptome,
);

open FF, "<$origTranscriptome";
open WW, ">$newTranscriptome";
while($line=<FF>){
	chomp($line);
	next if($line=~/^#/);
	@field = split(/\t/, $line);

	@attr = ();
	$field[8] = substr($field[8], 0, length($field[8]) - 1) if(substr($field[8], length($field[8]) - 1, 1) eq ";");
	@attr = split(/; /, $field[8]);
	$attrString = "";
	foreach $attr(@attr){
		next if($attr=~/reference_id/ or $attr=~/ref_gene_id/ or $attr=~/ref_gene_name/);
		$attrString .= $attr . "; ";
	}
	$field[8] = substr($attrString, 0, length($attrString) - 1);
	print WW join("\t", @field) . "\n";
}
close FF;
close WW;
