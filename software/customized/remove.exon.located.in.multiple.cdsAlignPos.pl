#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--inputExonPosInCodonAlign \\\n" .
                "--outputExonPosInSingleCodonAlign \\\n" .
                "--outputExonPosInMultipleCodonAlign \n";
	exit;
}

my ($inputExonPosInCodonAlign, $outputExonPosInSingleCodonAlign, $outputExonPosInMultipleCodonAlign);

GetOptions(
        'inputExonPosInCodonAlign=s'=>\$inputExonPosInCodonAlign,
        'outputExonPosInSingleCodonAlign=s'=>\$outputExonPosInSingleCodonAlign,
        'outputExonPosInMultipleCodonAlign=s'=>\$outputExonPosInMultipleCodonAlign,
);

my ($chr, $start, $stop, $strand, $cdsAlignPosList, @cdsAlignPos, $cdsAlignPos, %cdsAlignPosHash);
open SINGLE, ">$outputExonPosInSingleCodonAlign";
open MULTIPLE, ">$outputExonPosInMultipleCodonAlign";
open FF, "<$inputExonPosInCodonAlign";
# 1       23519   24451   +       AT1G01040:OG0011502:1-1284,AT1G01040:OG0011502:1-1284
while(my $line=<FF>){
#	print $line;
#	<STDIN>;
	chomp($line);
	($chr, $start, $stop, $strand, $cdsAlignPosList) =split(/\t/, $line);
	@cdsAlignPos = ();
	@cdsAlignPos = split(/,/, $cdsAlignPosList);
	%cdsAlignPosHash = ();
	foreach $cdsAlignPos(@cdsAlignPos){
		$cdsAlignPosHash{$cdsAlignPos} = 1;
	}
	@cdsAlignPos = ();
	@cdsAlignPos = keys(%cdsAlignPosHash);
	if($#cdsAlignPos > 0 ){
		print MULTIPLE $line . "\n";
	}else{
		print SINGLE join("\t", $chr, $start, $stop, $strand, $cdsAlignPos[0]) . "\n";
	}
}
close FF;
close MULTIPLE;
close SINGLE;
