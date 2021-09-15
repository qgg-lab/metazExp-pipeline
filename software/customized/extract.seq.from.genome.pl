#!/usr/bin/perl
use strict;
if($#ARGV<0){
	print $0 . " genome.fa\n";
	exit;
}
my ($genomeFasta) = $ARGV[0];
my (%genomeSeq, $line, $flag, @tt, $id, $seq);
my ($chr, $strand, $seqBegin, $seqEnd);
open FF, "<$genomeFasta";
while($line=<FF>){
	chomp($line);
	if($line=~/>(.*)/){
		$id = $1;
		@tt = ();
		@tt = split(/ / , $id);
		$id = $tt[0]
	}else{
		$genomeSeq{$id}.=$line;
	}
}
close FF;

$flag = "Y";
while(uc($flag) ne "N"){
	print "input \"chr,strand,begin,end:\"\n";
	$line=<STDIN>;
	chomp($line);
	($chr, $strand, $seqBegin, $seqEnd) = split(/,/, $line);
	$seq = substr($genomeSeq{$chr}, $seqBegin-1, $seqEnd - $seqBegin + 1);
	if($strand eq "-"){
		$seq=reverse($seq);
		$seq=~tr/ACGT/TGCA/;
	}
	print "seqLen: " . length("$seq") . "\n";
	print "seq:\n";
	print $seq . "\n";
	print "extract next seq in DNA?Y/N\n";
	$flag = <STDIN>;
	chomp($flag);
}
