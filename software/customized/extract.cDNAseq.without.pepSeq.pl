#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--inputCdnaSeqFile ./final.complete.trspt.cDNA.fa \\\n".
		"--inputPepSeqFile ./final.complete.trspt.pep.fa \\\n" . 
		"--outputCdnaWithoutPepSeqFile final.complete.trspt.cDNA.without.pep.fa \n";
	exit;
}
my ($outputDir);
my ($interproscan, $inputPepSeqFile, $inputCdnaSeqFile, $outputCdnaWithoutPepSeqFile);

GetOptions(
	'inputCdnaSeqFile=s'=>\$inputCdnaSeqFile,
	'inputPepSeqFile=s'=>\$inputPepSeqFile,
	'outputCdnaWithoutPepSeqFile=s'=>\$outputCdnaWithoutPepSeqFile,
);


my (%cDNAseq, %pepSeq, $id);
open FF, "<$inputCdnaSeqFile";
while(my $line =<FF>){
        chomp($line);
        if($line=~/>(.*) transcript_name:(.*) gene_id:(.*) gene_name:(.*)/){
                $id = $1;
        }else{
                $cDNAseq{$id}.=$line;
        }
}
close FF;

open FF, "<$inputPepSeqFile";
while(my $line =<FF>){
        chomp($line);
        if($line=~/>(.*) transcript_name:(.*) gene_id:(.*) gene_name:(.*) protein_id:(.*)/){
                $id = $1;
        }else{
                $pepSeq{$id}.=$line;
        }
}
close FF;

my @cDNAid=keys(%cDNAseq);
open WW, ">$outputCdnaWithoutPepSeqFile";
foreach $id(@cDNAid){
	if(not(exists($pepSeq{$id}))){
		print WW ">$id\n";
		print WW $cDNAseq{$id} . "\n";
	}
}
close WW;


