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
                "--psiDir ../013-gather-readCount-for-asAndAs/expts/AS/jcec\\\n" .
                "--controlExptIdListList SRX507922,SRX507923,SRX864604,SRX864605,SRX864609_SRX507922,SRX507923,SRX864604,SRX864605,SRX864609\\\n" .
                "--treatmentExptIdListList SRX335915,SRX335916,SRX335917,SRX335918,SRX335919_SRX5880876,SRX5880877,SRX5880878,SRX5891278,SRX5891279\\\n" .
		"--changedPsi 0.05\\\n" .
		"--diffQvalue 0.05 \\\n" .
		"--outputDir ./output.0001\\\n".
		"--finalOutputList 12345.39947.leaf.root.diff.splicing.as.by.psi.tsv_12345.39947.leaf.seed.diff.splicing.as.by.psi.tsv\n";
	exit;
}

my ($psiDir, $controlExptIdListList, $treatmentExptIdListList, $changedPsi, $diffQvalue, $finalOutputList, $outputDir, $minReadCoverage);
my ($jobId, $taxonId);
$jobId= time();

GetOptions(
	'jobId=s'=>\$jobId,
	'taxonId=s'=>\$taxonId,
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



#=== 以下是for(my $i=0;$i<=$#controlExptIdList;$i++){ =的循环体部分 ===
# 开启第i个比较任务


$controlExptIdList = $controlExptIdList[$i];
$treatmentExptIdList = $treatmentExptIdList[$i];
$finalOutput = $outputDir . "/" . $finalOutput[$i];

# 定义hash存放as的psi数据
my ($exptId, $asId, $psi, $countFile);
my (%asExptPsi, $asExptPsiHref, %asId, @asId);
%asId = ();
@asId = ();
%asExptPsi = ();
$asExptPsiHref=\%asExptPsi;

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

# 将expt中psi数读入asExptTpm
my (@field, $groupName, $psiFile, $asId, $IJC_SAMPLE, $SJC_SAMPLE, $IncFormLen, $SkipFormLen, $experimentId);
foreach $exptId(@controlExptId, @treatmentExptId){
	# 获得当前experiment对应的group
	$groupName = $exptToGroup{$exptId};
	# 当前experiment下基因表达psi
	$psiFile = $psiDir . "/" . $exptId;
	open FF, "<$psiFile";
	<FF>;
	# ASID    		IJC_SAMPLE_1    SJC_SAMPLE_1    IncFormLen      SkipFormLen     Experiment
	# OSJGA3SS0000000001    0       	6       	104     	100     	SRX335964
	# OSJGA3SS0000000002    36      	21      	105     	100     	SRX335964
	while(my $line=<FF>){	
		chomp($line);
		($asId, $IJC_SAMPLE, $SJC_SAMPLE, $IncFormLen, $SkipFormLen, $experimentId) = split(/\t/, $line);

		# 检测当前AS上在这个expt上覆盖的read数量是否达标，如果不达标，那么放弃该psi
		next if($IJC_SAMPLE + $SJC_SAMPLE < $minReadCoverage);

		# 计算当前AS在这个expt上的psi值
		$psi = ($IJC_SAMPLE/$IncFormLen)/($IJC_SAMPLE/$IncFormLen + $SJC_SAMPLE/$SkipFormLen);
		$asExptPsiHref->{$asId}->{$groupName} .= $psi . ",";
		$asId{$asId} = 1;
	}
	close FF;
}

# 生成as的psi列表矩阵
my $rndId=md5_hex($jobId . $taxonId . $controlExptIdList . $treatmentExptIdList);
my ($as_psi_matrix, $outputLine, @psi1, @psi2, $psiList1, $psiList2, %psiCL1, %psiCL2);
@asId = keys(%asId);
$as_psi_matrix = $outputDir . "/$jobId.$taxonId.$rndId.as.psi.matrix.tsv";
open WW, ">$as_psi_matrix";
print WW join("\t", "as_Id", "control", "treatment") . "\n";
# as_id              control              treatment
# Os12g0127400       4.5,6.3,6,2,2.3,3.1  9.3,12.3,3.4
foreach $asId(@asId){
	$outputLine = $asId;

	# 检测control中的psi值数量是否小于2，小于2则放弃该asId
	$psiList1 = substr($asExptPsiHref->{$asId}->{"control"}, 0, length($asExptPsiHref->{$asId}->{"control"}) - 1);
	@psi1 = ();
	@psi1 = split(/,/, $psiList1);
	next if($#psi1 + 1 < 2);
	# 检测control是否为常量
	%psiCL1 = ();
	foreach $psi(@psi1){
		$psiCL1{$psi}=1;
	}
	@psi1 = ();
	@psi1 = keys(%psiCL1);

	# 检测treatment中的psi值数量是否小于2，小于2则放弃该asId
	$psiList2 = substr($asExptPsiHref->{$asId}->{"treatment"}, 0, length($asExptPsiHref->{$asId}->{"treatment"}) - 1);
	@psi2 = ();
	@psi2 = split(/,/, $psiList2);
	next if($#psi2 + 1 < 2);
	# 检测treatment是否为常量
	%psiCL2 = ();
	foreach $psi(@psi2){
		$psiCL2{$psi}=1;
	}
	@psi2 = ();
	@psi2 = keys(%psiCL2);

	# 如果两个group中的psi值都是常量那么放弃该asId
	next if($#psi1 == 0 and $#psi2 == 0);

	# 输出该asId
	$outputLine .= "\t" . $psiList1 . "\t" . $psiList2;
	print WW $outputLine . "\n";
}
close WW;



# 调用R
my $R=Statistics::R->new(shared => 1);
my $cmd=<<EOF;
df<-read.table(file="$as_psi_matrix", sep="\\t", quote ="", header=T)
df\$changedPsi <- apply(df, 1, function(x) mean(as.numeric(unlist(strsplit(x[3], ","))))-mean(as.numeric(unlist(strsplit(x[2], ",")))))
df\$pvalue <- p.adjust(apply(df, 1, function(x) t.test(as.numeric(unlist(strsplit(x[2], ","))),as.numeric(unlist(strsplit(x[3], ","))), alternative='two.sided', paired = F)\$p.value), method="fdr")
rlt <- df[which((df[, "changedPsi"]>=$changedPsi | df[, "changedPsi"]<=-1 * $changedPsi) & df[, "pvalue"]< $diffQvalue), ]
write.table(rlt, file='$finalOutput', sep='\\t', quote=F, row.names = F)
EOF

$R->run($cmd);

$R->stop();

#print $cmd;
system("rm -rf $as_psi_matrix");
system("touch $outputDir/$jobId.$taxonId.$rndId.finished");

#=== 以上是for(my $i=0;$i<=$#controlExptIdList;$i++){ =的循环体部分 ===

} # for($i=0; $i<=$#controlExptIdList; $i++), 遍历所有比较
