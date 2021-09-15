#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--statisticOfSummarizedRlt \\\n" .
                "--classLablList \\\n" .
                "--taxonIdListInClass \\\n" .
		"--innerCdsExonWithConversedLevel \n";
	exit;
}

my ($exonPosInCodonAlign, @taxonId, $taxonId, $listOfSummarizedRlt, $statisticOfSummarizedRlt, 
	$classLablList, $taxonIdListInClass, $orthoCdsAlignPosConversedLevel);

GetOptions(
        'statisticOfSummarizedRlt=s'=>\$statisticOfSummarizedRlt,
	'classLablList=s'=>\$classLablList,
	'taxonIdListInClass=s'=>\$taxonIdListInClass,
	'orthoCdsAlignPosConversedLevel=s'=>\$orthoCdsAlignPosConversedLevel,
);

# 将物种谱系对应的taxonId读入
my ($i, $j);
my (%lineage, $lineageHref, @lineageLabel, $lineageLabel);
my (@taxonIdList);
$lineageHref=\%lineage;
@lineageLabel = split(/,/, $classLablList);
@taxonIdList = split(/#/, $taxonIdListInClass);
for($i=0; $i<=$#lineageLabel; $i++){
	$lineageHref->{$lineageLabel[$i]} = $taxonIdList[$i];
}

my (%orthoInnerCdsExonConserv, $orthoInnerCdsExonConservHref);
$orthoInnerCdsExonConservHref = \%orthoInnerCdsExonConserv;
my (@titleField, $titleField, @valueField, $valueField, %orthoCodonAlignPosToExonNumHash, $orthoCodonAlignPosToExonNumHashHref, $line);
$orthoCodonAlignPosToExonNumHashHref = \%orthoCodonAlignPosToExonNumHash;
my (@orthoCodonAlignPos, $orthoCodonAlignPos);
# 将exonStatistic读入，根据此orthoInnerCdsExon在各个物种中出现的次数
# 对其赋予保守水平
# orthoCodonAlignPos      3702    4577    39947   3847    4081
# OG0003417:879-1013      1       2       3       1       1
open FF, "<$statisticOfSummarizedRlt";
$line = <FF>;
chomp($line);
@titleField = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	# 将在每个物种中出现的exon数量存入orthoCodonAlignPosToExonNumHash
	@valueField = split(/\t/, $line);
	$orthoCodonAlignPos = $valueField[0];
	for($i=1; $i<=$#valueField; $i++){
		$orthoCodonAlignPosToExonNumHashHref->{$titleField[$i]} = $valueField[$i];
	}

	# 当前orthoCodonAlignPos的保守性级别对应的值
	# 列出所有谱系的标签label
	foreach $lineageLabel(@lineageLabel){
		# 提取对应谱系中所有的taxonId编号，存放在数组@taxonId中
		@taxonId = ();
		@taxonId = split(/,/, $lineageHref->{$lineageLabel});
		# 检查当前%orthoCodonAlignPosToExonNumHash中所有@taxonId是否都大于0
		my $checkRlt = 1;
		foreach $taxonId(@taxonId){
			$checkRlt *= $orthoCodonAlignPosToExonNumHashHref->{$taxonId};
		}
		if($checkRlt > 0){
			$orthoInnerCdsExonConservHref->{$orthoCodonAlignPos}->{$lineageLabel} = 1;
		}else{
			$orthoInnerCdsExonConservHref->{$orthoCodonAlignPos}->{$lineageLabel} = 0;
		}
	}
	
	# 当前orthoCodonAlignPos在各个物种内部是否保守，如果在2个以上位置出现则为保守，设置为1，否则设置为0
	foreach $titleField(@titleField){
		if($orthoCodonAlignPosToExonNumHashHref->{$titleField} >= 1){
			$orthoInnerCdsExonConservHref->{$orthoCodonAlignPos}->{$titleField} = $orthoCodonAlignPosToExonNumHashHref->{$titleField};
		}else{
			$orthoInnerCdsExonConservHref->{$orthoCodonAlignPos}->{$titleField} = 0;
		}
	}
}
close FF;

###########
# orthoCodonAlignPos在各个物种中的保守情况输出来
my @orthoCodonAlignPos = keys(%orthoInnerCdsExonConserv);
# 输出表头
my ($outputLine);
shift(@titleField);
open WW, ">$orthoCdsAlignPosConversedLevel";
print WW join("\t", "orthoCodonAlignPos", @lineageLabel, @titleField) . "\n";
foreach $orthoCodonAlignPos(@orthoCodonAlignPos){
	$outputLine = $orthoCodonAlignPos;
	foreach $lineageLabel(@lineageLabel){
		$outputLine .= "\t" . $orthoInnerCdsExonConservHref->{$orthoCodonAlignPos}->{$lineageLabel};
	}
	foreach $titleField(@titleField){
		$outputLine .= "\t" . $orthoInnerCdsExonConservHref->{$orthoCodonAlignPos}->{$titleField};
	}
	print WW $outputLine . "\n";
}
close WW;
