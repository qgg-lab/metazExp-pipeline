#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--inputExonPosInCdsSeq \\\n" .
		"--species \\\n" .
		"--outputExonPosInCdsSeq  \n";
	exit;
}

my ($inputExonPosInCdsSeq, $species, $outputExonPosInCdsSeq );

GetOptions(
        'inputExonPosInCdsSeq=s'=>\$inputExonPosInCdsSeq,
	'species=s'=>\$species,
        'outputExonPosInCdsSeq=s'=>\$outputExonPosInCdsSeq,
);

my ($line, @field, @posInCds, $posInCds, %geneId, @geneId);
# 1       3996    4276    +       Atha_AT1G01010_pep002:155-435,Atha_AT1G01010_pep001:155-435
# 1       4486    4605    +       Atha_AT1G01010_pep002:436-555
open FF, "<$inputExonPosInCdsSeq";
open WW, ">$outputExonPosInCdsSeq";
while($line=<FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);
	%geneId = ();
	@posInCds = ();
	@posInCds = split(/,/, $field[4]);
	foreach $posInCds(@posInCds){
		if($posInCds=~/($species)_(.*?)_pep\d+:\d+\-\d+/){
			$geneId{$2} = 1;		}
	}
	@geneId = ();
	@geneId = keys(%geneId);
#	print "@geneId";
#	<STDIN>;
	if($#geneId == 0){
		print WW $line . "\n";
	}
}
close FF;
close WW;
