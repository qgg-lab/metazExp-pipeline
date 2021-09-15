#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--species \"Bos taurus\" \\\n" .
		"--inputASfile  A5SS.catalog \\\n" . 
		"--inputAStype  A5SS \\\n" . 
		"--outputAsTsvFile 89462.A5ss.tsv \n\n";
	exit;
}

my ($species);
my ($inputASfile, $inputAStype, $outputAsTsvFile);

GetOptions(
        'species=s'=>\$species,
	'inputASfile=s'=>\$inputASfile, 
	'inputAStype=s'=>\$inputAStype, 
	'outputAsTsvFile=s'=>\$outputAsTsvFile,
);

my ($line, @fields, $tableFields, $values, $sql);
my ($asId, $asType, $chr, $strand, $start, $end, $firstExonEnd, $firstExonStart_0base, $secondExonEnd, $secondExonStart_0base, $downstreamEE, $downstreamES, $exonEnd, $exonStart_0base, $flankingEE, $flankingES, $longExonEnd, $longExonStart_0base, $riExonEnd, $riExonStart_0base, $shortEE, $shortES, $upstreamEE, $upstreamES, $geneId, $geneSymbol, $firstAltExonSeries, $firstIsoformIdList, $firstIsoformNum, $secondAltExonSeries, $secondIsoformIdList, $secondIsoformNum, $pfam, $go, $variantPointNum, $variantPointTypeCmb, $asOrthId, $conservedSpeciesNum, $jcecExperimentNum, $jcExperimentNum, $discoveryApproach);

$asType = $inputAStype;

open WW, ">$outputAsTsvFile";
open FF, "<$inputASfile";
	$line = <FF>;
	chomp($line);
	@fields = ();
	@fields = split(/\t/, $line);

while($line=<FF>){
	chomp($line);
	@fields = ();
	@fields = split(/\t/, $line);
	
	$asType = $inputAStype;
	$species = $species;

	($start, $end, $firstExonEnd, $firstExonStart_0base, $secondExonEnd, $secondExonStart_0base, $downstreamEE, $downstreamES, $exonEnd, $exonStart_0base, $flankingEE, $flankingES, $longExonEnd, $longExonStart_0base, $riExonEnd, $riExonStart_0base, $shortEE, $shortES, $upstreamEE, $upstreamES, $firstIsoformNum, $secondIsoformNum, $variantPointNum, $conservedSpeciesNum, $jcecExperimentNum, $jcExperimentNum) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

	($chr, $strand, $asId, $geneId, $geneSymbol, $firstAltExonSeries, $firstIsoformIdList, $secondAltExonSeries, $secondIsoformIdList, $pfam, $go, $variantPointTypeCmb, $asOrthId, $discoveryApproach) = ("", "+", "", "", "", "", "", "", "", "", "", "", "", "");
	
	if(uc($asType) eq "A5SS"){

		#ASID    GeneID  geneSymbol      chr     strand  longExonStart_0base     longExonEnd     shortES shortEE flankingES      flankingEE
		($asId, $geneId, $geneSymbol, $chr, $strand, $longExonStart_0base, $longExonEnd, $shortES, $shortEE, $flankingES, $flankingEE) = ($fields[0], $fields[1], $fields[2], $fields[3], $fields[4], $fields[5], $fields[6], $fields[7], $fields[8], $fields[9], $fields[10]);

	}elsif(uc($asType) eq "A3SS"){

		#ASID    GeneID  geneSymbol      chr     strand  longExonStart_0base     longExonEnd     shortES shortEE flankingES      flankingEE
		($asId, $geneId, $geneSymbol, $chr, $strand, $longExonStart_0base, $longExonEnd, $shortES, $shortEE, $flankingES, $flankingEE) = ($fields[0], $fields[1], $fields[2], $fields[3], $fields[4], $fields[5], $fields[6], $fields[7], $fields[8], $fields[9], $fields[10]);

        }elsif(uc($asType) eq "SE"){
		
		#ASID    GeneID  geneSymbol      chr     strand  exonStart_0base exonEnd upstreamES      upstreamEE      downstreamES    downstreamEE
		($asId, $geneId, $geneSymbol, $chr, $strand, $exonStart_0base, $exonEnd, $upstreamES, $upstreamEE, $downstreamES, $downstreamEE) = ($fields[0], $fields[1], $fields[2], $fields[3], $fields[4], $fields[5], $fields[6], $fields[7], $fields[8], $fields[9], $fields[10]);

	}elsif(uc($asType) eq "RI"){
		
		#ASID    GeneID  geneSymbol      chr     strand  riExonStart_0base       riExonEnd       upstreamES      upstreamEE      downstreamES    downstreamEE
		($asId, $geneId, $geneSymbol, $chr, $strand, $riExonStart_0base, $riExonEnd, $upstreamES, $upstreamEE, $downstreamES, $downstreamEE) = ($fields[0], $fields[1], $fields[2], $fields[3], $fields[4], $fields[5], $fields[6], $fields[7], $fields[8], $fields[9], $fields[10]);

	}elsif(uc($asType) eq "MXE"){
		
		#ASID GeneID geneSymbol chr strand 1stExonStart_0base 1stExonEnd 2ndExonStart_0base 2ndExonEnd upstreamES upstreamEE downstreamES downstreamEE
		($asId, $geneId, $geneSymbol, $chr, $strand, $firstExonStart_0base, $firstExonEnd, $secondExonStart_0base, $secondExonEnd, $upstreamES, $upstreamEE, $downstreamES, $downstreamEE) = ($fields[0], $fields[1], $fields[2], $fields[3], $fields[4], $fields[5], $fields[6], $fields[7], $fields[8], $fields[9], $fields[10], $fields[11], $fields[12]);
	}

	$geneId=~s/"//g;
	$geneSymbol=~s/"//g;
	if($chr=~/^chr(.*)/){
		$chr=$1;
	}

	print WW "asId, asType, species, chr, strand, start, end, 1stExonEnd, 1stExonStart_0base, 2ndExonEnd, 2ndExonStart_0base, downstreamEE, downstreamES, exonEnd, exonStart_0base, flankingEE, flankingES, longExonEnd, longExonStart_0base, riExonEnd, riExonStart_0base, shortEE, shortES, upstreamEE, upstreamES, geneID, geneSymbol, firstAltExonSeries, firstIsoformIdList, firstIsoformNum, secondAltExonSeries, secondIsoformIdList, secondIsoformNum, pfam, go, variantPointNum, variantPointTypeCmb, asOrthId, conservedSpeciesNum, jcecExperimentNum, jcExperimentNum, discoveryApproach" . "_____" . "\"$asId\", \"$asType\", \"$species\", \"$chr\", \"$strand\", $start, $end, $firstExonEnd, $firstExonStart_0base, $secondExonEnd, $secondExonStart_0base, $downstreamEE, $downstreamES, $exonEnd, $exonStart_0base, $flankingEE, $flankingES, $longExonEnd, $longExonStart_0base, $riExonEnd, $riExonStart_0base, $shortEE, $shortES, $upstreamEE, $upstreamES, $geneId, $geneSymbol, \"$firstAltExonSeries\", \"$firstIsoformIdList\", $firstIsoformNum, \"$secondAltExonSeries\", \"$secondIsoformIdList\", $secondIsoformNum, \"$pfam\", \"$go\", $variantPointNum, \"$variantPointTypeCmb\", \"$asOrthId\", $conservedSpeciesNum, $jcecExperimentNum, $jcExperimentNum, \"$discoveryApproach\"\n";
}

close FF;
close WW;
