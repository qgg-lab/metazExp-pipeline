#!/usr/bin/perl
use strict;
use Digest::MD5 qw(md5 md5_hex md5_base64);
use Getopt::Long;
use Statistics::R;
use DBI;
if($#ARGV<0){
my $textContent=<<EOF;
比如有2组比较：
第1组：controlExptIdList=ERX2120647,ERX2120648,ERX2120649
       treatmentExptIdList=ERX2120650,ERX2120651,ERX2120652
       finalOutput=12345.39947.root.leaf.diff.expr.as.by.DESeq2.tsv
第2组：controlExptIdList=ERX2120647,ERX2120648,ERX2120649
       treatmentExptIdList=ERX2120649,ERX2120648,ERX2120650
       finalOutput=12345.39947.root.callus.diff.expr.as.by.DESeq2.tsv

这里需要把2组比较合并在一起，用";"将2个comparison的control/treatmentExptIdList分隔开来
	controlExptIdListList=controlExptIdList1;controlExptIdList2
	treatmentExptIdListList=treatmentExptIdList1;treatmentExptIdList2
	finalOutputList=finalOutput1;finalOutput2
具体参数设置如下
EOF
	print $textContent;
	print "\nperl $0 \\\n" . 
		"--taxonId 39947\\\n" .
		"--jobId 12345 \\\n" .
		"--pytho /mnt/home/liujind1/software/Python-2.7.13/bin/python \\\n" .
		"--rMatStat /mnt/home/liujind1/software/rMATS-STAT/MATS_LRT.py \\\n" .
		"--detectPackageFromSingleCompare /mnt/home/liujind1/software/customized/detect.diff.splicing.as.by.rmatsStat.pl \\\n".
                "--psiDir ../013-gather-readCount-for-asAndAs/expts/AS/jcec\\\n" .
                "--controlExptIdListList SRX507922,SRX507923,SRX864604,SRX864605,SRX864609_SRX507922,SRX507923,SRX864604,SRX864605,SRX864609\\\n" .
                "--treatmentExptIdListList SRX335915,SRX335916,SRX335917,SRX335918,SRX335919_SRX5880876,SRX5880877,SRX5880878,SRX5891278,SRX5891279\\\n" .
		"--changedPsi 0.05 \\\n" .
		"--minReadCoverage 5 \\\n" .
		"--diffQvalue 0.05 \\\n" .
		"--outputDir ./output.0001\\\n".
		"--finalOutputList 12345.39947.leaf.root.diff.splicing.as.by.psi.tsv_12345.39947.leaf.seed.diff.splicing.as.by.psi.tsv\n";
	exit;
}

my ($psiDir, $controlExptIdListList, $treatmentExptIdListList, $changedPsi);
my ($diffQvalue, $finalOutputList, $outputDir, $minReadCoverage, $python, $rMatStat);
my ($jobId, $taxonId);
$jobId= time();

GetOptions(
	'jobId=s'=>\$jobId,
	'taxonId=s'=>\$taxonId,
	'python=s'=>\$python,
	'rMatStat=s'=>\$rMatStat,
        'psiDir=s'=>\$psiDir,
        'controlExptIdListList=s'=>\$controlExptIdListList,
        'treatmentExptIdListList=s'=>\$treatmentExptIdListList,
        'changedPsi=s'=>\$changedPsi,
	'minReadCoverage=s'=>\$minReadCoverage,
        'diffQvalue=s'=>\$diffQvalue,
	'outputDir=s'=>\$outputDir,
        'finalOutputList=s'=>\$finalOutputList,
);

# 创建输出目录
system("mkdir -p $outputDir");

# 生成比较数组，包括输入的control和treatment的ExptId的列表 和 对应的输出文件名
my (@controlExptIdList, $controlExptIdList, @treatmentExptIdList, $treatmentExptIdList, @finalOutput, $finalOutput);
@controlExptIdList = split(/_/, $controlExptIdListList);
@treatmentExptIdList = split(/_/, $treatmentExptIdListList);
@finalOutput = split(/_/, $finalOutputList);

# 循环列出每一组比较的输入和输出
for(my $i=0;$i<=$#controlExptIdList;$i++){

	$controlExptIdList = $controlExptIdList[$i];
	$treatmentExptIdList = $treatmentExptIdList[$i];
	$finalOutput = $outputDir . "/" . $finalOutput[$i];

	# 调用rMatStat
	my $cmd= "perl $detectPackageFromSingleCompare " . 
		" --python $python " . 
		" --rmatsStat $rMatStat" . 
		" --asReadCountDir $psiDir" . 
		" --controlExptIdList $controlExptIdList" . 
		" --treatmentExptIdList $treatmentExptIdList" .
		" --diff_cutoff $changedPsi" .
		" --threadNum $threadNum" . 
		" --outputDir $finalOutput";

	system("mv $finalOutput $finalOutput.backup");
	system("mv $finalOutput.backup/rMATS_Result_P.txt $finalOutput");
	system("rm -rf $finalOutput.backup");
	system("touch $outputDir/$jobId.$taxonId.$rndId.finished");


} # for($i=0; $i<=$#controlExptIdList; $i++), 遍历所有比较
