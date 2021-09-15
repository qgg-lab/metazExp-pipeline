#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--inputAsWithOrthExonFile A5SS.catalog.with.orthTag \\\n" .
		"--mappedCoordFile mapped.A5SS.coords.in.hg38.tsv \\\n" .
		"--outputAsFileInHg38 A5SS.cataloag.with.coords.in.hg38 \n";
	exit;
}

my ($inputAsWithOrthExonFile, $mappedCoordFile, $outputAsFileInHg38);
GetOptions(
	'inputAsWithOrthExonFile=s'=>\$inputAsWithOrthExonFile, 
	'mappedCoordFile=s'=>\$mappedCoordFile, 
	'outputAsFileInHg38=s'=>\$outputAsFileInHg38, 
);

my ($line, @fields, $field);

# read AS
my (@fieldTitle, $i,  %as, $orthId, $asCoord, $asId, $coordOrthExonId, $coordName);
open FF, "<$mappedCoordFile";
# chr11   77891839        77891840        GGALA5SS0000011669:ORTH0000078644:flankingES    0       -
while($line=<FF>){
	chomp($line);
	@fields = ();
	@fields = split(/\t/, $line);
	($asId, $coordOrthExonId, $coordName) = split(/:/, $fields[3]);
	${$as{$asId . "#" . $coordName}}{"hg38Coord"} = join("#", $fields[0], $fields[1], $fields[2], $fields[5]);
	${$as{$asId . "#" . $coordName}}{"orthExonId"} = $coordOrthExonId;
}

my ($asIdCoordName);
open FF, "<$inputAsWithOrthExonFile";
$line=<FF>;
#ASID GeneID geneSymbol chr strand longExonStart_0base longExonEnd shortES shortEE flankingES flankingEE
chomp($line);
@fieldTitle = ();
@fieldTitle = split(/\t/, $line);
open WW, ">$outputAsFileInHg38";
print WW $line . "\n";
# GGALA5SS0000018079 "ENSGALG00000000164" "MYBPH"    chr26 -  chr11:77891839:77891840 chr11:70988698:70988699 - - - -
while($line=<FF>){
	#GGALA5SS0000018079 "ENSGALG00000000164" "MYBPH"    chr26 -      ORTH0000196814 ORTH0000196814 - - - -
	chomp($line);
	@fields = ();
	@fields = split(/\t/, $line);
	print WW join("\t", $fields[0], $fields[1], $fields[2], $fields[3], $fields[4]);
	for($i=5; $i<=$#fields; $i++){
		if($fields[$i]=~/ORTH/){
			$asIdCoordName = $fields[0] . "#" . $fieldTitle[$i];
			if(exists(${$as{$asIdCoordName}}{"hg38Coord"}) and ${$as{$asIdCoordName}}{"orthExonId"} eq $fields[$i]){
				print WW "\t" . ${$as{$asIdCoordName}}{"hg38Coord"};
			}else{
				print WW "\t-";
			}
		}else{
			print WW "\t" . $fields[$i];
		}
	}
	print WW "\n";
}
close WW;
close FF;
