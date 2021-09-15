#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--blastRlt \\\n" .
                "--outputBestHit \n";
	exit;
}

my ($blastRlt, $outputBestHit);

GetOptions(
        'blastRlt=s'=>\$blastRlt,
        'outputBestHit=s'=>\$outputBestHit,
);

my (%bestHit, $href, $line, @field, $geneId, $hitGeneId, $evalue);
$href = \%bestHit;
open FF, "<$blastRlt";
# LOC108855588    AT2G31470       46.94   392     190     9       5       392     10      387     2e-108    327
# LOC108805682    AT3G19510       63.11   740     198     19      3       688     4       722     0.0       785
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	($geneId, $hitGeneId, $evalue) = ($field[0], $field[1], $field[10]);
	if(not exists($bestHit{$geneId})){
		$href->{$geneId}->{"hitGeneId"} = $hitGeneId;
		$href->{$geneId}->{"evalue"} = $evalue;
	}elsif($evalue<$href->{$geneId}->{"evalue"}){
		$href->{$geneId}->{"hitGeneId"} = $hitGeneId;
		$href->{$geneId}->{"evalue"} = $evalue;
	}
}
close FF;

my @geneId = keys(%bestHit);
open WW, ">$outputBestHit";
foreach $geneId(@geneId){
	print WW join("\t", $geneId, $href->{$geneId}->{"hitGeneId"}) . "\n";
}
close WW;


