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
       finalOutput=12345.39947.root.leaf.diff.expr.gene.by.DESeq2.tsv
第2组：controlExptIdList=ERX2120647,ERX2120648,ERX2120649
       treatmentExptIdList=ERX2120649,ERX2120648,ERX2120650
       finalOutput=12345.39947.root.callus.diff.expr.gene.by.DESeq2.tsv

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
                "--tpmDir ../013-gather-readCount-for-geneAndAs/expts/gene/\\\n" .
                "--controlExptIdListList SRX507922,SRX507923,SRX864604,SRX864605,SRX864609_SRX507922,SRX507923,SRX864604,SRX864605,SRX864609\\\n" .
                "--treatmentExptIdListList SRX335915,SRX335916,SRX335917,SRX335918,SRX335919_SRX5880876,SRX5880877,SRX5880878,SRX5891278,SRX5891279\\\n" .
		"--diffFoldchange 2\\\n" .
		"--diffQvalue 0.05 \\\n" .
		"--outputDir ./mult.tpm\\\n".
		"--finalOutputList 12345.39947.leaf.root.diff.expr.gene.by.DESeq2.tsv_12345.39947.leaf.seed.diff.expr.gene.by.DESeq2.tsv\n";
	exit;
}

my ($tpmDir, $controlExptIdListList, $treatmentExptIdListList, $diffFoldchange, $diffQvalue, $finalOutputList, $outputDir);
my ($jobId, $taxonId);
$jobId= time();

GetOptions(
	'jobId=s'=>\$jobId,
	'taxonId=s'=>\$taxonId,
        'tpmDir=s'=>\$tpmDir,
        'controlExptIdListList=s'=>\$controlExptIdListList,
        'treatmentExptIdListList=s'=>\$treatmentExptIdListList,
        'diffFoldchange=s'=>\$diffFoldchange,
        'diffQvalue=s'=>\$diffQvalue,
	'outputDir=s'=>\$outputDir,
        'finalOutputList=s'=>\$finalOutputList,
);

# 创建输出目录
system("mkdir -p $outputDir");

# 换算foldchange为log2
$diffFoldchange = log($diffFoldchange)/log(2);

# 生成比较数组，包括输入的control和treatment的ExptId的列表 和 对应的输出文件名
my (@controlExptIdList, $controlExptIdList, @treatmentExptIdList, $treatmentExptIdList, @finalOutput, $finalOutput);
@controlExptIdList = split(/_/, $controlExptIdListList);
@treatmentExptIdList = split(/_/, $treatmentExptIdListList);
@finalOutput = split(/_/, $finalOutputList);

# 循环列出每一组比较的输入和输出
for(my $i=0;$i<=$#controlExptIdList;$i++){



#=== 以下是for(my $i=0;$i<=$#controlExptIdList;$i++){ =的循环体部分 ===
# 开启第i个比较任务


$controlExptIdList = $controlExptIdList[$i];
$treatmentExptIdList = $treatmentExptIdList[$i];
$finalOutput = $outputDir . "/" . $finalOutput[$i];

# 定义hash存放gene的tpm数据
my ($exptId, $geneId, $tpm, $countFile);
my (%geneExptTpm, $geneExptTpmHref, %geneId, @geneId);
%geneId = ();
@geneId = ();
%geneExptTpm = ();
$geneExptTpmHref=\%geneExptTpm;

# 拆解获得control和treatment中每个experiment的编号
my (@controlExptId, @treatmentExptId, $exptId);
my (%exptToGroup);
@controlExptId = ();
@controlExptId = split(/,/, $controlExptIdList);
@treatmentExptId = ();
@treatmentExptId = split(/,/, $treatmentExptIdList);
# 建立从experiment到control和treatment关系
foreach $exptId(@controlExptId){
	$exptToGroup{$exptId} = "control";
}
foreach $exptId(@treatmentExptId){
	$exptToGroup{$exptId} = "treatment";
}

# 将expt中tpm数读入geneExptTpm
my (@field, $groupName, $tpmFile);
foreach $exptId(@controlExptId, @treatmentExptId){
	# 获得当前experiment对应的group
	$groupName = $exptToGroup{$exptId};
	# 当前experiment下基因表达tpm
	$tpmFile = $tpmDir . "/" . $exptId;
	open FF, "<$tpmFile";
	<FF>;
	# geneId  readCount       Coverage        FPKM    TPM
	# Os12g0624100    0       0.0     0.0     0.0
	# Os03g0243750    0       0.000000        0.000000        0.000000
	# Os11g0263000    429     9.699306        2.590241        4.017003
	while(my $line=<FF>){	
		chomp($line);
		@field = split(/\t/, $line);
		$geneId = $field[0];
		$tpm = $field[4];
		$geneExptTpmHref->{$geneId}->{$groupName} .= $tpm . ",";
		$geneId{$geneId} = 1;
	}
	close FF;
}

# 生成gene的tpm列表矩阵
my $rndId=md5_hex($jobId . $taxonId . $controlExptIdList . $treatmentExptIdList);
my ($gene_tpm_matrix, $outputLine);
@geneId = keys(%geneId);
$gene_tpm_matrix = $outputDir . "/$jobId.$taxonId.$rndId.gene.tpm.matrix.tsv";
open WW, ">$gene_tpm_matrix";
print WW join("\t", "gene_Id", "control", "treatment") . "\n";
# gene_id            control              treatment
# Os12g0127400       4.5,6.3,6,2,2.3,3.1  9.3,12.3,3.4
foreach $geneId(@geneId){
	$outputLine = $geneId;
	$outputLine .= "\t" . substr($geneExptTpmHref->{$geneId}->{"control"}, 0, length($geneExptTpmHref->{$geneId}->{"control"}) - 1);
	$outputLine .= "\t" . substr($geneExptTpmHref->{$geneId}->{"treatment"}, 0, length($geneExptTpmHref->{$geneId}->{"treatment"}) - 1);
	print WW $outputLine . "\n";
}
close WW;

# 调用DESeq2
my $R=Statistics::R->new(shared => 1);
my $cmd=<<EOF;
df<-read.table(file="$gene_tpm_matrix", sep="\\t", quote ="", header=T)
df\$foldChange <- apply(df, 1, function(x) mean(as.numeric(unlist(strsplit(x[3], ","))))/mean(as.numeric(unlist(strsplit(x[2], ",")))))
df\$pvalue <- p.adjust(apply(df, 1, function(x) t.test(as.numeric(unlist(strsplit(x[2], ","))),as.numeric(unlist(strsplit(x[3], ","))), alternative='two.sided', paired = F)\$p.value), method="fdr")
rlt <- df[which((df[, "foldChange"]>=$diffFoldchange | df[, "foldChange"]<=1/$diffFoldchange) & df[, "pvalue"]< $diffQvalue), ]
write.table(rlt, file='$finalOutput', sep='\\t', quote=F, row.names = F)
EOF

$R->run($cmd);

$R->stop();

#print $cmd;
system("rm -rf $gene_tpm_matrix");
system("touch $outputDir/$jobId.$taxonId.$rndId.finished");

#=== 以上是for(my $i=0;$i<=$#controlExptIdList;$i++){ =的循环体部分 ===

} # for($i=0; $i<=$#controlExptIdList; $i++), 遍历所有比较
