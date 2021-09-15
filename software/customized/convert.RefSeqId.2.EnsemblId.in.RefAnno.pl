#!/usr/bin/perl
use strict;
use Getopt::Long;
if($#ARGV < 0){
	print $0 . " \\\n" . 
		"\t\t --refAnnoFile       ref_ARS-UCD1.2_top_level.gff3      \\\n" .
		"\t\t --seqIdMapFile      GenomeSeqIdMapping.tsv             \\\n" .
		"\t\t --outputNewAnnoFile    new.ref_ARS-UCD1.2_top_level.gff3  \\\n";
	exit(0);
}

my (%refIdMapEnsemblId);
my ($refAnnoFile, $seqIdMapFile, $outputNewAnnoFile);
GetOptions(
        'refAnnoFile=s'=>\$refAnnoFile,
        'seqIdMapFile=s'=>\$seqIdMapFile,
        'outputNewAnnoFile=s'=>\$outputNewAnnoFile,
);

my (@tt, $line);

open FF, "<$seqIdMapFile";
#NW_020192282.1  NKLS02002198.1
#NW_020191789.1  NKLS02001705.1
while($line = <FF>){
	chomp($line);
	@tt = ();
	@tt = split(/\t/, $line);
	$refIdMapEnsemblId{$tt[0]} = $tt[1];
}
close FF;

open FF, "<$refAnnoFile";
open WW, ">$outputNewAnnoFile";

while($line = <FF>){
	chomp($line);
	@tt = ();
	@tt = split(/\t/, $line);
	if($#tt == 8 and exists($refIdMapEnsemblId{$tt[0]}) and $refIdMapEnsemblId{$tt[0]} ne "-"){
		print WW join("\t", $refIdMapEnsemblId{$tt[0]}, $tt[1], $tt[2], $tt[3], $tt[4], $tt[5], $tt[6], $tt[7], $tt[8]) . "\n";
	}
}
close FF;
close WW;
