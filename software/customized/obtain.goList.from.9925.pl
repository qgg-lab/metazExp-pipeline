#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--origGOlist \\\n" .
                "--blastRlt \\\n" .
		"--output \n";
	exit;
}

my ($origGOlist, $queryProteomeFasta, $blastRlt, $output);
my (%origGeneToGoList, $geneId, $goList, %queryGeneToOrigGene, $line, $evalue, @tmp);
my $qToOHref = \%queryGeneToOrigGene;

GetOptions(
        'origGOlist=s'=>\$origGOlist,
        'blastRlt=s'=>\$blastRlt,
        'output=s'=>\$output,
);

open FF, "<$origGOlist";
while($line=<FF>){
	chomp($line);
	($geneId, $goList)  = split(/\t/, $line);
	$origGeneToGoList{$geneId} = $goList;
}
close FF;

open FF, "<$blastRlt";
# ENSOARG00020013517      ENSCHIG00000011913      96.129  310     12      0       1       310     1       310     0.0     602
# ENSOARG00020000404      ENSCHIG00000013739      96.353  521     18      1       1       520     56      576     0.0     1011
while($line=<FF>){
	chomp($line);
	@tmp = split(/\t/, $line);
	if(not exists($qToOHref->{$tmp[0]})){
		$qToOHref->{$tmp[0]}->{"geneId"}= $tmp[1];
		$qToOHref->{$tmp[0]}->{"evalue"}= $tmp[10];
	}elsif($qToOHref->{$tmp[0]}->{"evalue"}>= $tmp[10]){
		$qToOHref->{$tmp[0]}->{"geneId"}= $tmp[1];
		$qToOHref->{$tmp[0]}->{"evalue"}= $tmp[10];
	}
}
close FF;

my ($targetGeneId);
open WW, ">$output";
my @geneId = keys(%queryGeneToOrigGene);
foreach $geneId(@geneId){
	$targetGeneId = $qToOHref->{$geneId}->{"geneId"};
	if(exists($origGeneToGoList{$targetGeneId})){
		print WW join("\t", $geneId, $origGeneToGoList{$targetGeneId}) . "\n";
	}
}
close WW;
