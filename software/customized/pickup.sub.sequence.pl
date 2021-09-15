#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--genomeSeqFileList 9031.genome.fa,9796.genome.fa,9823.genome.fa,9913.genome.fa,9940.genome.fa \\\n" .
                "--genomeLabelList 9031,9796,9823,9913,9940\\\n" .
		"--exonFileList 9031.exon.bed,9796.exon.bed,9823.exon.bed,9913.exon.bed,9940.exon.bed \\\n" .
		"--orthExonListFile final.orthExon.matrix.tsv \n";
	exit;
}

my ($genomeSeqFileList, @genomeSeqFile, $genomeSeqFile);
my ($genomeLabelList, @genomeLabel, $genomeLabel);
my ($exonFileList, @exonFile, $exonFile);
my ($orthExonListFile);
GetOptions(
        'genomeSeqFileList=s'=>\$genomeSeqFileList,
        'genomeLabelList=s'=>\$genomeLabelList,
	'exonFileList=s'=>\$exonFileList,
	'orthExonListFile=s'=>\$orthExonListFile,
);

my (%genomeSeq, $i, $line, @tt, $id);
my (%exon);
# read genome sequence into hash
@genomeSeqFile = split(/,/, $genomeSeqFileList);
@genomeLabel = split(/,/, $genomeLabelList);
@exonFile = split(/,/, $exonFileList);
for($i=0; $i<=$#genomeSeqFile; $i++){
	
	$genomeLabel = $genomeLabel[$i];
	$genomeSeqFile = $genomeSeqFile[$i];
	$exonFile = $exonFile[$i];

	open FF, "<$genomeSeqFile";
	while($line = <FF>){
        	chomp($line);
	        if($line=~/>/){
                	@tt = ();
        	        @tt = split(/ /, $line);
	                if($tt[0]=~/>ref\|(.*?)\|/ or $tt[0]=~/>gi.*\|ref\|(.*?)\|/ or  $tt[0]=~/>(.*)/){
                	        $id = $1;
        	        }
	        }else{
                	${$genomeSeq{$genomeLabel}}{$id} .=uc($line);
        	}
	}
	close FF;
	print "Finish load " . $genomeLabel . " genome.\n";

	open FF, "<$exonFile";
	while($line=<FF>){
		# chr1    5272    5524    Exon1   0       -
		chomp($line);
		@tt = ();
		@tt = split(/\t/, $line);
		${$exon{$genomeLabel}}{$tt[3]} = $line;
	}
	close FF;
	print "Finish load " . $genomeLabel . " exons.\n";;
}

my (%orthExon, @fieldTitle, @field, $orthExonId);
open FF, "<$orthExonListFile";
$line=<FF>;
chomp($line);
@fieldTitle = ();
@fieldTitle = split(/\t/, $line);
shift(@fieldTitle);

while($line=<FF>){
	#orthId  9031    9796    9823    9913    9940
	#ORTH0000000002  Exon116411      Exon225128      Exon118201      Exon213018      Exon185323
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);
	$orthExonId = shift(@field);
	for(my $i=0; $i<=$#field; $i++){
		${$orthExon{$orthExonId}}{$fieldTitle[$i]} = $field[$i];
#		print join("\t", $orthExonId, ${$orthExon{$orthExonId}}{$fieldTitle[$i]}) . "\n";
#		<STDIN>;
	}
}
close FF;
print "Finish load orthExonMatrix\n";


my ($inputFlag, $region, $chr, $chain, $start, $stop, $st, $score, $exonId);
$inputFlag = "Y";
while($inputFlag eq "Y"){
	while(1){
		print "Input orthExonId:\n";
		$orthExonId = <STDIN>;
		chomp($orthExonId);
		if(exists($orthExon{$orthExonId})){
			last;
		}else{
			print $orthExonId . " does not exists!\n";
		}
	}
	
	foreach $genomeLabel(@genomeLabel){
		print join("\t", $genomeLabel,${$exon{$genomeLabel}}{${$orthExon{$orthExonId}}{$genomeLabel}}) . "\n";
	}

	foreach $genomeLabel(@genomeLabel){
		($chr, $start, $stop, $exonId, $score, $chain) = ("", "", "", "", "", "");
		($chr, $start, $stop, $exonId, $score, $chain) = split(/\t/, ${$exon{$genomeLabel}}{${$orthExon{$orthExonId}}{$genomeLabel}});
		if($chr=~/chr0*(.*)/){
			$chr = $1;
		}
		$st = substr(${$genomeSeq{$genomeLabel}}{$chr}, $start-1, $stop - $start +1);
		if($chain eq "-"){
			$st = reverse($st);
			$st=~tr/ACGT/TGCA/;
		}
		print $genomeLabel . ":" . $st . "\n";
	}

	print "Do you want next orthExon sequences?\n";
	$inputFlag = <STDIN>;
	chomp($inputFlag);
	$inputFlag = uc($inputFlag);	
}
