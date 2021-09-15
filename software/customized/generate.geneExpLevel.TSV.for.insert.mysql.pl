#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--inputGeneExpListFile origin.geneExp.list.tsv \\\n" .
		"--species \"Gallus gallus\" \\\n" .
                "--taxon 9031 \\\n" .
                "--outputGeneExpTsvFile geneExp.tsv\n";
	exit;
}

my ($inputGeneExpListFile, $taxon, $outputGeneExpTsvFile, $species);

GetOptions(
        'inputGeneExpListFile=s'=>\$inputGeneExpListFile,
        'taxon=s'=>\$taxon,
	'species=s'=>\$species,
        'outputGeneExpTsvFile=s'=>\$outputGeneExpTsvFile,
);

my ($line, $i, @fieldTitle, $fieldTitle, @field, $field);
open WW, ">$outputGeneExpTsvFile";
open FF, "<$inputGeneExpListFile";
# geneID  geneName        Chromo  Strand  Start   End     Coverage        FPKM    TPM     ExperimentId
# ENSGALG00000045335      -       W       -       2988277 2989653 1.899027        0.142906        0.415449        SRX1036607
<FF>;
while($line =<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	print WW "species, taxon, geneId, geneName, chromo, strand, start, end, coverage, fpkm, tpm, experimentId" . "_____";
	print WW "\"" . $species . "\", \"" . $taxon . "\", \"" . $field[0] . "\", \"" . $field[1] . "\", \"" . $field[2] . "\", \"" . $field[3] . "\", " . $field[4] . ", " . $field[5] . ", " . $field[6] . ", " . $field[7] . ", " . $field[8] . ", \"" . $field[9] . "\"\n";
}
close FF;
