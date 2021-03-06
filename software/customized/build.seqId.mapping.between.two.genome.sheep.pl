#!/usr/bin/perl
use strict;
use Getopt::Long; 
if($#ARGV < 0){
	print $0 . " \\\n" .
		"\t\t --firstGenomeFile  ensemblGenome.fa \\\n" . 
		"\t\t --secondGenomeFile refSeqGenome.fa  \\\n" .
		"\t\t --outputSeqIdMappingFile GenomeSeqIdMapping.tsv \\\n";
	exit(0);
}
my ($firstGenomeFile, $secondGenomeFile, $outputSeqIdMappingFile);
my (%firstGenomeSeq, %secondGenomeSeq);

GetOptions(
        'firstGenomeFile=s'=>\$firstGenomeFile,
        'secondGenomeFile=s'=>\$secondGenomeFile,
        'outputSeqIdMappingFile=s'=>\$outputSeqIdMappingFile,
);

my ($line, $id, @tt);

open FF, "<$firstGenomeFile";
#>NKLS02001232.1 dna_sm:primary_assembly primary_assembly:ARS-UCD1.2:NKLS02001232.1:1:1531:1 REF
while($line = <FF>){
	chomp($line);
	if($line=~/>/){
		@tt = ();
		@tt = split(/ /, $line);
		$id = substr($tt[0], 1);
	}else{
		$firstGenomeSeq{$id} .=uc($line);
	}
}
close FF;

open FF, "<$secondGenomeFile";
#>gi|417531923|ref|NC_019467.1| Ovis
#>gi|417531923|ref|NC_019467.1| Ovis aries breed Texel chromosome 10, Oar_v3.1, whole genome shotgun sequence
while($line = <FF>){
	chomp($line);
	if($line=~/>/){
		@tt = ();
		@tt = split(/ /, $line);
		if($tt[0]=~/>gi\|\d+\|ref\|(.*?)\|/){
			$id = $1;
		#	print $id . "\n";
		}
	}else{
		$secondGenomeSeq{$id} .=uc($line);
	}
}
close FF;

my ($flag);
open WW, ">$outputSeqIdMappingFile";
my @firstGenomeSeqId = keys %firstGenomeSeq;
my @secondGenomeSeqId = keys %secondGenomeSeq;
for(my $i=0; $i<=$#secondGenomeSeqId; $i++){
	$flag = 0;
	for(my $j=0; $j<=$#firstGenomeSeqId; $j++){
		if($secondGenomeSeq{$secondGenomeSeqId[$i]} eq $firstGenomeSeq{$firstGenomeSeqId[$j]}){
			print WW join("\t", $secondGenomeSeqId[$i], $firstGenomeSeqId[$j]) . "\n";
			$flag = 1;
			last;
		}
	}
	if($flag==0){
		print WW join("\t", $secondGenomeSeqId[$i], "-") . "\n";
	}
}
close WW;
