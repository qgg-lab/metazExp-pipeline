#!/usr/bin/perl
use strict;
my ($genomeFasta, $newGenomeFasta) = @ARGV;
if($#ARGV<1){
	print "$0 genome.fa new.genome.fa\n";
	exit;
}

# 读入genome序列
my (%genomeSeq, $genomeSeqHref);
$genomeSeqHref=\%genomeSeq;
my ($line, @chr, $chr);
my $timeLabel = localtime();
open FF, "<$genomeFasta";
while($line=<FF>){
        chomp($line);
        if($line=~/>(.*)/){
                @chr = split(/ /, $1);
                $chr = $chr[0];
        }else{
                $genomeSeqHref->{$chr}.=$line;
        }
}
close FF;

open WW, ">$newGenomeFasta";
@chr = keys(%genomeSeq);
foreach $chr(@chr){
	print WW ">$chr\n";
	print WW $genomeSeq{$chr} . "\n";
}
close WW;
