#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--exonPosInCodonAlign \\\n" .
		"--taxonId \\\n" .
		"--lineageList \\\n" .
                "--orthoCdsAlignPosConversedLevel \\\n" .
                "--outputInnerCdsExonWithConservLevel \n";
	exit;
}

my ($exonPosInCodonAlign, $taxonId, $lineageList, $orthoCdsAlignPosConversedLevel, $outputInnerCdsExonWithConservLevel);

GetOptions(
        'exonPosInCodonAlign=s'=>\$exonPosInCodonAlign,
	'taxonId=s'=>\$taxonId,
	'lineageList=s'=>\$lineageList,
        'orthoCdsAlignPosConversedLevel=s'=>\$orthoCdsAlignPosConversedLevel,
        'outputInnerCdsExonWithConservLevel=s'=>\$outputInnerCdsExonWithConservLevel,
);

my ($line, $i, $j, @lineage, $lineage);
my (%orthoCdsAlignPos, $orthoCdsAlignPosHref, @titleField, @valueField);
my (@orthoCdsAlignPos, $orthoCdsAlignPos);
$orthoCdsAlignPosHref = \%orthoCdsAlignPos;
@lineage=split(/_/, $lineageList);
# 将orthoCdsAlignPos读入hash 
open FF, "<$orthoCdsAlignPosConversedLevel";
# orthoCodonAlignPos      Angi    Mono    Dico    Poal    Capp    Sola    Faba    3702    4577    39947   3847    4081
# OG0007697:754-879       1       1       1       1       1       1       1       0       0       0       1       1
# OG0003417:879-1013      1       1       1       1       1       1       1       0       1       1       0       0 
$line = <FF>;
chomp($line);
@titleField = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@valueField = split(/\t/, $line);
	$orthoCdsAlignPos = $valueField[0];
	for($i=1; $i<=$#valueField; $i++){
		$orthoCdsAlignPosHref->{$orthoCdsAlignPos}->{$titleField[$i]} = $valueField[$i];
	}
	$orthoCdsAlignPosHref->{$orthoCdsAlignPos}->{"conservLevel"} = "-";
	foreach $lineage(@lineage){
		if($orthoCdsAlignPosHref->{$orthoCdsAlignPos}->{$lineage} == 1){
			$orthoCdsAlignPosHref->{$orthoCdsAlignPos}->{"conservLevel"} = $lineage;
			last;
		}
	}
}
close FF;

# 读取物种的innerCdsExon，重新调整一下
# geneId  exonId  chr     exonStart       exonStop        exonStrand      orthoCdsAlignPos
# AT1G01020       AthaExon000001  1       7157    7232    -       orthoPos108404
# AT1G01020       AthaExon000002  1       7384    7450    -       orthoPos124639

# 调整后
# geneId	exonId		chr	start	stop	strand	orthoCdsAlignPos	conservLevel	duplication
# AT1G01020	AthaExon000001  1       7157    7232    -       orthoPos108404	        Angi/Mono/Dico	2
my ($exonId, $chr, $start, $stop, $strand, $geneId, $orthoCdsAlignPos, $conserveLevel);
open FF, "<$exonPosInCodonAlign";
<FF>;
open WW, ">$outputInnerCdsExonWithConservLevel";
print WW join("\t", "geneId", "exonId", "chr", "start", "stop", "strand", "orthoCdsAlignPos", "conservLevel", "duplication") . "\n";
while($line=<FF>){
	chomp($line);
	($geneId, $exonId, $chr, $start, $stop, $strand, $orthoCdsAlignPos) = split(/\t/, $line);
	print WW join("\t", $geneId, $exonId, $chr, $start, $stop, $strand, $orthoCdsAlignPos, $orthoCdsAlignPosHref->{$orthoCdsAlignPos}->{"conservLevel"}, $orthoCdsAlignPosHref->{$orthoCdsAlignPos}->{$taxonId}) . "\n";
}
close WW;
