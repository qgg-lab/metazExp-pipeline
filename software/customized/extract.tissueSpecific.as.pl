#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--taxonId \\\n" .
		"--asTableMysqlTsv \n";
	exit;
}

my ($asTableMysqlTsv, $taxonId);

GetOptions(
	'taxonId=s'=>\$taxonId,
	'asTableMysqlTsv=s'=>\$asTableMysqlTsv,
);

my ($specialHighNum, $specialLowNum) = (0, 0);
my ($line, $i, $valueList, $nameList, @value, @name, $value, $name, %as);
open FF, "<$asTableMysqlTsv";
# asId, asType, chr, strand, start, stop, orthAsGroupId, asAnnoSource, geneId, geneSymbol, genePfamList, geneGoList, otherAsIdListInSameGene, 1stExonEnd, 1stExonStart_0base, 2ndExonEnd, 2ndExonStart_0base, downstreamEE, downstreamES, exonEnd, exonStart_0base, flankingEE, flankingES, longExonEnd, longExonStart_0base, riExonEnd, riExonStart_0base, shortEE, shortES, upstreamEE, upstreamES, inclExonSeries, inclTrsptIdList, inclTrsptIdListInOrigAnnoGtf, inclTrsptIdListInImprAnnoGtf, exclExonSeries, exclTrsptIdList, exclTrsptIdListInOrigAnnoGtf, excluTrsptIdListInImprAnnoGtf, editSiteInBothAltSeq, cutSizeInInclAlt, insertSizeInInclAlt, cutSizeInExclAlt, insertSizeInExclAlt, psiSpecialHighTissueName, psiSpecialLowTissueName, mirnaTargetList, validedByEnoughExpts, 1stExon, 2ndExon, upstreamExon, downstreamExon, flankingExon, exon, longExon, shortExon, riExon, inclusionAltSeq, exclusionAltSeq___SLYCMXE0000004044, MXE, 1, -, 92879885, 92882550, NA, RNAseqMapping, Solyc01g104530.3, NA, PF00069, GO:0004672,GO:0005524,GO:0006468, SLYCA3SS0000004294,SLYCSE0000010969,SLYCRI0000001697,SLYCSE0000017082,SLYCSE0000021950,SLYCSE0000042202, 92881041, 92880981, 92881324, 92881251, 92882550, 92882392, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 92879995, 92879884, 92882393..92882550,92880982..92881041,92879885..92879995, NA, NA, NA, 92882393..92882550,92881252..92881324,92879885..92879995, NA, NA, NA, 158, 60, 73, 73, 60, NA, NA, NA, No, GCAACGCAAATGAATAACATCAAGTCTTGCAAAGGAACTGCATTCTGGATGGCCCCAGAG, GGACATTAAGTGTGCCAATATATTGGTGGATGCTAATGGTTCAGTGAAGCTTGCAGATTTTGGACTGGCAAAG, GTTGTTAATAGGAAGAGCAATGGATATGGTACTCCTGCTGATATATGGAGTTTGGGTTGCACTGTCTTAGAGATGTTAACTGGACAAATTCCTTATTCTCACTTGGAAGGG, GACGAGAGCAAACTGTATATCTTTCTTGAGCTTGTTACAAAAGGTTCACTTGCTAGTGTCTATCGCAAGTATCGCTTGCGTGATTCCCATGTTTCAGATTACACCAGGCAAATCTTGAGTGGGTTGCATTATCTTCATTCCAGAGAAGTGATGCACAG, NA, NA, NA, NA, NA, GACGAGAGCAAACTGTATATCTTTCTTGAGCTTGTTACAAAAGGTTCACTTGCTAGTGTCTATCGCAAGTATCGCTTGCGTGATTCCCATGTTTCAGATTACACCAGGCAAATCTTGAGTGGGTTGCATTATCTTCATTCCAGAGAAGTGATGCACAGGCAACGCAAATGAATAACATCAAGTCTTGCAAAGGAACTGCATTCTGGATGGCCCCAGAGGTTGTTAATAGGAAGAGCAATGGATATGGTACTCCTGCTGATATATGGAGTTTGGGTTGCACTGTCTTAGAGATGTTAACTGGACAAATTCCTTATTCTCACTTGGAAGGG, GACGAGAGCAAACTGTATATCTTTCTTGAGCTTGTTACAAAAGGTTCACTTGCTAGTGTCTATCGCAAGTATCGCTTGCGTGATTCCCATGTTTCAGATTACACCAGGCAAATCTTGAGTGGGTTGCATTATCTTCATTCCAGAGAAGTGATGCACAGGGACATTAAGTGTGCCAATATATTGGTGGATGCTAATGGTTCAGTGAAGCTTGCAGATTTTGGACTGGCAAAGGTTGTTAATAGGAAGAGCAATGGATATGGTACTCCTGCTGATATATGGAGTTTGGGTTGCACTGTCTTAGAGATGTTAACTGGACAAATTCCTTATTCTCACTTGGAAGGG
while($line=<FF>){
	chomp($line);
	($nameList, $valueList) = split(/___/, $line);
	#print $valueList;
	#<STDIN>;
	#print $nameList;
	#<STDIN>;
	@value = split(/, /, $valueList);
	@name = split(/, /, $nameList);
	for($i=0; $i<=$#name; $i++){
	#	print "name:".$name[$i];
	#	<STDIN>;
	#	print "value:" . $value[$i];
	#	<STDIN>;
		$as{$name[$i]} = $value[$i];
	}
	#print join("\t", "psiSpecialHighTissueName:" . $as{"psiSpecialHighTissueName"}, "psiSpecialLowTissueName:". $as{"psiSpecialLowTissueName"});
	#<STDIN>;
	$specialHighNum++ if($as{"psiSpecialHighTissueName"} ne "NA");
	$specialLowNum++ if($as{"psiSpecialLowTissueName"} ne "NA");
}
close FF;

print join("\t", $taxonId, $specialHighNum, $specialLowNum) . "\n";
