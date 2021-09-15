#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--trsptTagFile trspt.with.orign.tissue.exptList.tsv\\\n" .
                "--filterFile filter.trspt.by.tissue.conf \\\n" .
                "--outputRmTrsptIdListFile removed.trsptId.list.1 \n";
	exit;
}

my ($trsptTagFile, $filterFile, $outputRmTrsptIdListFile);

GetOptions(
        'trsptTagFile=s'=>\$trsptTagFile,
        'filterFile=s'=>\$filterFile,
        'outputRmTrsptIdListFile=s'=>\$outputRmTrsptIdListFile,
);

my ($line, @tmp, $i);
# 将过滤标准文件读入到hash中
my (%tissueToExptNum);
open FF, "<$filterFile";
# tissue  tissueCutoff
# ear     2
# embryo  3
# total       210
<FF>;
while($line=<FF>){
	chomp($line);
	@tmp = ();
	@tmp = split(/\t/, $line);	
	$tissueToExptNum{$tmp[0]}=$tmp[1];
}
close FF;

# 读取每个转录本，只关注新组装的转录本，查看其是否符合条件。如果不符合条件，那么输出该转录本的编号
my ($trsptId, $tissue, $exptList, @exptId, $delFlag, $tissueId,  $tissueExptList, $totalExptNum, @tt);
open WW, ">$outputRmTrsptIdListFile";
open FF, "<$trsptTagFile";
# SRX901422.27396.6       StringTie       shoot:SRX901422
# Zm00001d013026_T001     gramene seedling:SRX5794451,SRX5794453,SRX5794454       stem:SRX5387719,SRX5387724
while($line=<FF>){

	chomp($line);
	@tmp = split(/\t/, $line);
	$trsptId = $tmp[0];

	next if($tmp[1] ne "StringTie");

	# 默认删除，符合条件时将delFlag修改为0，表示不删除
	$delFlag = 1;
	$totalExptNum = 0;

	for($tissueId=2; $tissueId<=$#tmp; $tissueId++){
		$tissueExptList = $tmp[$tissueId];
		@tt = ();
		@tt = split(/:/, $tissueExptList);
		$tissue = $tt[0];
		$exptList = $tt[1];
		@exptId = ();
		@exptId = split(/,/, $exptList);
		if(exists($tissueToExptNum{$tissue}) and $#exptId+1 >= $tissueToExptNum{$tissue}){
			$delFlag = 0;
			last;
		}else{
			$totalExptNum+=$#exptId+1;
		}
	}	

	# 检查是否存在符合的条件
	if($delFlag == 1 and $totalExptNum < $tissueToExptNum{"total"}){
		print WW $trsptId . "\n";
	}
}
close FF;
close WW;
