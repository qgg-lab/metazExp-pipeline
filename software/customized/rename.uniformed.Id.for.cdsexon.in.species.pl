#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--exonPosInCodonAlign \\\n" .
                "--abbr \\\n" .
                "--uniformedExonList \n";
	exit;
}

my ($exonPosInCodonAlign, $uniformedExonList, $abbr);

GetOptions(
        'exonPosInCodonAlign=s'=>\$exonPosInCodonAlign,
        'abbr=s'=>\$abbr,
        'uniformedExonList=s'=>\$uniformedExonList,
);

# 读取exonPos in codon alignment
my ($exonNum);
open WW, ">$uniformedExonList";
open FF, "<$exonPosInCodonAlign";
# 1       7157    7232    -       AT1G01020:OG0006459:676-802
# 1       7384    7450    -       AT1G01020:OG0006459:609-675
# 1       7564    7649    -       AT1G01020:OG0006459:514-608
while(my $line=<FF>){
	$exonNum++;
	print WW $abbr . "Exon" . sprintf("%06d", $exonNum) . "\t" . $line;
}
close FF;
close WW;
