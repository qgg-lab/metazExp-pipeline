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
                "--countDir ../013-gather-readCount-for-geneAndAs/expts/gene/\\\n" .
                "--controlExptIdListList ERX2120647,ERX2120648,ERX2120649_ERX2120647,ERX2120648,ERX2120649\\\n" .
                "--treatmentExptIdListList ERX2120650,ERX2120651,ERX2120652_ERX2120649,ERX2120648,ERX2120650\\\n" .
		"--diffFoldchange 2\\\n" .
		"--diffQvalue 0.05 \\\n" .
		"--outputDir ./\\\n".
		"--finalOutputList 12345.39947.root.leaf.diff.expr.gene.by.DESeq2.tsv_12345.39947.root.callus.diff.expr.gene.by.DESeq2.tsv\n";
	exit;
}

my ($countDir, $controlExptIdListList, $treatmentExptIdListList, $diffFoldchange, $diffQvalue, $finalOutputList, $outputDir);
my ($jobId, $taxonId);
$jobId= time();

GetOptions(
	'jobId=s'=>\$jobId,
	'taxonId=s'=>\$taxonId,
        'countDir=s'=>\$countDir,
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

# 定义hash存放gene的readcount数据
my ($exptId, $geneId, $count, $countFile);
my (%geneExptCount, $geneExptCountHref, %geneId, @geneId);
%geneId = ();
@geneId = ();
%geneExptCount = ();
$geneExptCountHref=\%geneExptCount;

# 拆解获得control和treatment中每个experiment的编号
my (@controlExptId, @treatmentExptId);
@controlExptId = ();
@controlExptId = split(/,/, $controlExptIdList);
@treatmentExptId = ();
@treatmentExptId = split(/,/, $treatmentExptIdList);

# 将expt中read数读入geneExptCount
my @field;
foreach $exptId(@controlExptId, @treatmentExptId){
	$countFile = $countDir . "/" . $exptId;
	open FF, "<$countFile";
	<FF>;
	# geneId  readCount       Coverage        FPKM    TPM
	# Os12g0624100    0       0.0     0.0     0.0
	# Os03g0243750    0       0.000000        0.000000        0.000000
	# Os11g0263000    429     9.699306        2.590241        4.017003
	while(my $line=<FF>){	
		chomp($line);
		@field = split(/\t/, $line);
		$geneId = $field[0];
		$count = $field[1];
		$geneExptCountHref->{$geneId}->{$exptId} = $count;
		$geneId{$geneId} = 1;
	}
	close FF;
}

# 生成gene的read count矩阵
my $rndId=md5_hex($jobId . $taxonId . $controlExptIdList . $treatmentExptIdList);
my ($gene_count_matrix, $outputLine);
@geneId = keys(%geneId);
$gene_count_matrix = $outputDir . "/$jobId.$taxonId.$rndId.gene.count.matrix.tsv";
open WW, ">$gene_count_matrix";
print WW join(",", "gene_Id", @controlExptId, @treatmentExptId) . "\n";
# gene_id,ERX2120647,ERX2120648,ERX2120649,ERX2120650,ERX2120651,ERX2120652
# Os12g0127400,45,63,68,28,23,31
foreach $geneId(@geneId){
	$outputLine = $geneId;
	foreach $exptId(@controlExptId, @treatmentExptId){
		if(exists($geneExptCountHref->{$geneId}->{$exptId})){
			$outputLine .= "," . $geneExptCountHref->{$geneId}->{$exptId};
		}else{
			$outputLine .= ",0";
		}
	}
	print WW $outputLine . "\n";
}
close WW;

# 生成样本信息表
my $sampleInfo = $outputDir . "/$jobId.$taxonId.$rndId.sample.info.tsv";
open WW, ">$sampleInfo";
#exptId  group
#ERX2120647      control
print WW join("\t", "exptId", "group") . "\n";
foreach $exptId(@controlExptId){
	print WW join("\t", $exptId, "control") . "\n";
}
foreach $exptId(@treatmentExptId){
	print WW join("\t", $exptId, "treatment") . "\n";
}
close WW;

# 调用DESeq2
my $R=Statistics::R->new(shared => 1);
my $cmd=<<EOF;
library(DESeq2)
countData <-as.matrix(read.csv("$gene_count_matrix", row.names="gene_Id"))
colData <- read.csv("$sampleInfo", sep="\\t", row.names=1)
countData <- countData[, rownames(colData)]
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~group)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res\$padj), ]
padjRlt <- resOrdered[which(resOrdered[, "padj"]<=$diffQvalue),]
filterRlt1 <- padjRlt[which(padjRlt[, "log2FoldChange"]>=$diffFoldchange),]
filterRlt2 <- padjRlt[which(padjRlt[, "log2FoldChange"]<=-$diffFoldchange),]
finalRlt <- rbind(filterRlt1, filterRlt2)
write.table(finalRlt, file="$finalOutput", sep="\\t", row.names=T, col.names=T)
EOF
$R->run($cmd);
$R->stop();
#print $cmd;
#system("rm -rf $gene_count_matrix $sampleInfo");
system("touch $outputDir/$jobId.$taxonId.$rndId.finished");

#=== 以上是for(my $i=0;$i<=$#controlExptIdList;$i++){ =的循环体部分 ===

} # for($i=0; $i<=$#controlExptIdList; $i++), 遍历所有比较
