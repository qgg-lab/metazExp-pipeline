#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--exonPosInCodonAlign \\\n" .
                "--checkRlt \n";
	exit;
}

my ($exonPosInCodonAlign, $checkRlt);

GetOptions(
        'exonPosInCodonAlign=s'=>\$exonPosInCodonAlign,
        'checkRlt=s'=>\$checkRlt,
);
# 读取mergedInnerCdsExon及其对应在orthoCdsAlignPos
my ($line, @field, @orthoCdsAlignPos, $orthoCdsAlignPos, %orthoCdsAlignPos, @tt);
open WW, ">$checkRlt";
open FF, "<$exonPosInCodonAlign";
# 1       7157    7232    -       AT1G01020:OG0006459:802-883,AT1G01020:OG0006459:802-883,AT1G01020:OG0006459:802-883,AT1G01020:OG0006459:802-883,AT1G01020:OG0006459:802-883,AT1G01020:OG0006459:802-883,AT1G01020:OG0006459:802-883,AT1G01020:OG0006459:802-883 
while($line=<FF>){
	chomp($line);
	@field = ();
	@field = split(/\t/, $line);
	%orthoCdsAlignPos = ();
	@orthoCdsAlignPos = ();
	@orthoCdsAlignPos = split(/,/, $field[4]);
	foreach $orthoCdsAlignPos(@orthoCdsAlignPos){
		$orthoCdsAlignPos{$orthoCdsAlignPos} = 1;
	}
	@tt = ();
	@tt = keys(%orthoCdsAlignPos);
	if($#tt==0){
		pop(@field);
		print WW join("\t", @field, $tt[0]) . "\n";
	}
}
close FF;
close WW;
