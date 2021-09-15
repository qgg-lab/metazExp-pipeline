#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--orthoGroupOfLongestPep \\\n" .
                "--uniqPepFastaList \\\n" .
		"--seqPathOfOrthoGroupWithFullUniqPep \n";
	exit;
}

my ($orthoGroupOfLongestPep, $uniqPepFastaList, $seqPathOfOrthoGroupWithFullUniqPep);

GetOptions(
        'orthoGroupOfLongestPep=s'=>\$orthoGroupOfLongestPep,
        'uniqPepFastaList=s'=>\$uniqPepFastaList,
        'seqPathOfOrthoGroupWithFullUniqPep=s'=>\$seqPathOfOrthoGroupWithFullUniqPep,
);
# 将uniqPepFasta中uniqPep序列读入到uniqPep哈希中，同时将uniqPepId聚类到geneIdToUniqPepIdList哈希中
my (%geneIdToUniqPepIdList, %uniqPep, $geneId, $uniqPepId, $line);
my (@uniqPepFasta, $uniqPepFasta);
@uniqPepFasta = split(/,/, $uniqPepFastaList);
foreach $uniqPepFasta(@uniqPepFasta){
	open FF, "<$uniqPepFasta";
	# >Atha_AT1G09170_pep004
	# >Gmax_GLYMA_06G240600_pep002
	while($line=<FF>){
		chomp($line);
		if($line=~/>(.*?)_(.*)_(pep\d+)/){
			$uniqPepId = substr($line, 1);
			$geneId = $2;
			if(not exists($geneIdToUniqPepIdList{$geneId})){
				$geneIdToUniqPepIdList{$geneId} = $uniqPepId;
			}else{
				$geneIdToUniqPepIdList{$geneId}.="," . $uniqPepId;
			}
			$uniqPep{$uniqPepId} = "";
		}else{
			$uniqPep{$uniqPepId} .= $line;
		}
	}
	close FF;
}

# 读取orthoGroup.tsv中的OG编号和OG编号对应的longest uniqPepId
# 解析出longest uniqPepId列表所对应的geneId列表
# 逐个提取geneId，然后在哈希geneIdToUniqPepIdList中找到geneId对应的uniqPepId，再用uniqPepId在%uniqPep中获得pep序列
my ($ogId, @uniqPepIdList, @uniqPepId, $uniqPepIdList, $uniqPepId, @geneId, %geneId);
open FF, "<$orthoGroupOfLongestPep";
<FF>;
# Orthogroup Atha                  Gmax                         Osat                     Slyc                         Zmay
# OG0000000  Atha_AT2G21200_pep001 Gmax_GLYMA_04G006300_pep001  Osat_Os08g0118500_pep001 Slyc_Solyc01g110630.3_pep001 Zmay_Zm00001d014682_pep00
while($line=<FF>){
	chomp($line);
	@uniqPepIdList = ();
	@uniqPepIdList = split(/\t/, $line);
	$ogId = shift(@uniqPepIdList);
	%geneId = ();
	@geneId = ();
	# 为该ogId创建fasta文件
	open WW, ">" . $seqPathOfOrthoGroupWithFullUniqPep . "/" . $ogId . ".fa";
	# 找出longestUniqPep对应的geneId，存放在@geneId中
	foreach $uniqPepIdList(@uniqPepIdList){
		@uniqPepId = ();
		@uniqPepId = split(/,/, $uniqPepIdList);
		foreach $uniqPepId(@uniqPepId){
			if($uniqPepId=~/(.*?)_(.*)_(pep\d+)/){
				$geneId = $2;
				$geneId{$geneId}=1;
			}
		}
	}
	@geneId = ();
	@geneId = keys(%geneId);
	
	# 找出geneId对应的所有的uniqPepId，然后输出uniqPepId对应的序列
	foreach $geneId(@geneId){
		@uniqPepId = ();
		@uniqPepId = split(/,/, $geneIdToUniqPepIdList{$geneId});
		foreach $uniqPepId(@uniqPepId){
			print WW ">" . $uniqPepId . "\n";
			print WW $uniqPep{$uniqPepId} . "\n";
		}
	}
	close WW;
}
close FF;
