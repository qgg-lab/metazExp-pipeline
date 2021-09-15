#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--uniquePepFasta \\\n" .
                "--trsptCdsStrucGtf \\\n" .
                "--consistencyCheckRltTsv \n";
	exit;
}

my ($uniquePepFasta, $trsptCdsStrucGtf, $consistencyCheckRltTsv);

GetOptions(
        'uniquePepFasta=s'=>\$uniquePepFasta,
        'trsptCdsStrucGtf=s'=>\$trsptCdsStrucGtf,
        'consistencyCheckRltTsv=s'=>\$consistencyCheckRltTsv,
);

my (%trsptCdsStruc, $trsptCdsStrucHref, %pepToTrsptIdList);
$trsptCdsStrucHref=\%trsptCdsStruc;
# 将trspt的cds的外显子串结构以字符串形式读入到hash中
my ($cmd, $exonLineText, @exonLine, $exonLine, @field, $geneId, $trsptId);
$cmd = "grep -P \"\\texon\\t\" $trsptCdsStrucGtf";
$exonLineText = `$cmd`;
@exonLine = split(/\n/, $exonLineText);
foreach $exonLine(@exonLine){
	@field = ();
	@field = split(/\t/, $exonLine);
	if($field[8]=~/geneId ".*?"; transcript_id "(.*?)"; /){
		$trsptId = $1;
	}else{
		if(not exists($trsptCdsStrucHref->{$trsptId})){
			$trsptCdsStrucHref->{$trsptId} = $field[3] . ".." . $field[4];
		}else{
			$trsptCdsStrucHref->{$trsptId} .= "," . $field[3] . ".." . $field[4];
		}
	}
}

# 将uniq pep fasta中序列名依次读取出来，比较内部trsptId对应exon是否相同
# >AT5G09410___SRX399568.18675.3XXXATHAA3SS0000000398#AT5G09410.2XXXATHASE0000000241#SRX1660582.15918.3XXXATHAA3SS0000000399#AT5G09410.1
my ($geneIdAndTrsptIdListText, @geneIdAndTrsptIdList, $geneIdAndTrsptIdList);
my ($tsptIdList, @trsptId, $trsptId, $geneId, %exonStru, @exonStru);
$cmd = "grep \">\" " . $uniquePepFasta . " |awk -F \'>\' \'{print \$2}\' ";
$geneIdAndTrsptIdListText = `$cmd`;
@geneIdAndTrsptIdList = split(/\n/, $geneIdAndTrsptIdListText);
open WW, ">$consistencyCheckRltTsv";
print WW join("\t", "pepId", "exonStrucNum") . "\n";
foreach $geneIdAndTrsptIdList(@geneIdAndTrsptIdList){
	($geneId, $tsptIdList) = split(/___/, $geneIdAndTrsptIdList);
	@trsptId = ();
	@trsptId = split(/#/, $tsptIdList);
	%exonStru = ();
	# 利用hash对exonStruc进行聚类
	foreach $trsptId(@trsptId){
		$exonStru{$trsptCdsStrucHref->{$trsptId}} += 1;
	}
	# 判断exonStruc种类数，从而确定一个pep对应的exonStruc是否来自一个
	@exonStru = ();
	@exonStru =keys(%exonStru);
	if($#exonStru>0){
		print WW join("\t", $geneIdAndTrsptIdList, $#exonStru+1) . "\n";
	}
}
close WW;
