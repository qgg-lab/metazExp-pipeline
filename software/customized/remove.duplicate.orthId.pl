#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--inputOrthMap mapped.gal6.exon.bed \\\n" .
                "--inputOriginOrthCoord orthExon.coord.tsv \\\n" .
                "--remainingOrthCoord 1.orthExon.coord.tsv\n";
	exit;
}

my ($inputOrthMap, $inputOriginOrthCoord, $remainingOrthCoord);

GetOptions(
        'inputOrthMap=s'=>\$inputOrthMap,
        'inputOriginOrthCoord=s'=>\$inputOriginOrthCoord,
        'remainingOrthCoord=s'=>\$remainingOrthCoord,
);

my (%orthMapping, $line, @fields,  $hitPosition, %keptedOrth, @hitPosition, $orthId);
open FF, "<$inputOrthMap";
#chr1    147889380       147889478       ORTH0000000026  0       +
#chr5    38147753        38147876        ORTH0000000033  0       -
#chr5    34403052        34403128        ORTH0000000034  0       +
while($line=<FF>){
	chomp($line);
	@fields = ();
	@fields = split(/\t/, $line);
	$hitPosition = join("#", $fields[0], $fields[1], $fields[2], $fields[5]);
	${$orthMapping{$hitPosition}}{"orthIdList"}.=$fields[3] . "#";
	${$orthMapping{$hitPosition}}{"orthIdNum"}++;
}
close FF;

@hitPosition = keys(%orthMapping);
foreach $hitPosition(@hitPosition){
	if(${$orthMapping{$hitPosition}}{"orthIdNum"} == 1){
		$orthId = ${$orthMapping{$hitPosition}}{"orthIdList"};
		$orthId = substr($orthId, 0, length($orthId)-1);
		$keptedOrth{$orthId} = 1;
	}
}

open FF, "<$inputOriginOrthCoord";
open WW, ">$remainingOrthCoord";
#chr2    164915111       164915273       ORTH0000000009  0       -
#chr9    123169901       123170046       ORTH0000000016  0       -
while($line=<FF>){
	chomp($line);
	@fields = ();
	@fields = split(/\t/, $line);
	if(exists($keptedOrth{$fields[3]})){
		print WW $line . "\n";
	}
}
close FF;
close WW;
