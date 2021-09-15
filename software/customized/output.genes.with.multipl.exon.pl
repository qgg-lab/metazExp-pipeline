#!/usr/bin/perl
use strict;
my ($geneTrsptId)=@ARGV;
#Zm00001d027230   Zm00001d027230_T001
#Zm00001d027230   Zm00001d027230_T001
#Zm00001d027230   Zm00001d027230_T001
# 每一行对一个外显子，第1列为基因编号，第2列为转录本编号

# 统计方法：
# trspt->gene的hash
# gene -> multipleExon/singleExon的hash
# 识别trspt是否为multipleExon，如果是那么通过trspt->gene修改gene的类型
my (%gene, %trspt, @trsptId, $trsptId, @geneId, $geneId);
my (@tt);
open FF, "<$geneTrsptId";
while(my $line=<FF>){
	chomp($line);
	@tt = split(/\t/, $line);	
	$gene{$tt[0]} = 0;
	${$trspt{$tt[1]}}{"gene"} = $tt[0];
	${$trspt{$tt[1]}}{"exonNum"}+=1;
}
close FF;

@trsptId = keys(%trspt);
foreach $trsptId(@trsptId){
	$gene{${$trspt{$trsptId}}{"gene"}} = 1 if(${$trspt{$trsptId}}{"exonNum"}>1);
}

@geneId = keys(%gene);
foreach $geneId(@geneId){
	print $geneId . "\n" if($gene{$geneId}==1);
}
