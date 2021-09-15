#!/usr/bin/perl
use strict;
use Getopt::Long;
use List::Util qw/sum/;
use Statistics::R;

if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--threadNum \\\n" .
		"--multiCmpProgram multi.compare.find.diff.expr.gene.by.DESeq2.pl\\\n" .
                "--taxonId 39947 \\\n" .
                "--jobId 1234 \\\n" .
                "--exptDir ../013/expts/gene \\\n" .
                "--exptGroupInfo sample.info.tsv \\\n" .
                "--statTest Yes \\\n" .
                "--foldChange 2\\\n" .
                "--pvalue 0.05 \\\n" .
                "--minReplicateNum 3\\\n" .
                "--minGroupNum 5 \\\n" .
		"--outputDir specificityByDeseq2Dir\n";
	exit;
}

my ($taxonId, $jobId, $exptDir, $exptGroupInfo, $threadNum, $multiCmpProgram,
$statTest, $foldChange, $pvalue, 
$minReplicateNum, $minGroupNum, $outputDir);

GetOptions(
	'threadNum=s'=>\$threadNum,
	'multiCmpProgram=s'=>\$multiCmpProgram,
        'taxonId=s'=>\$taxonId,
        'jobId=s'=>\$jobId,
        'exptDir=s'=>\$exptDir,
        'exptGroupInfo=s'=>\$exptGroupInfo,
	'statTest=s'=>\$statTest,
	'foldChange=s'=>\$foldChange,
	'pvalue=s'=>\$pvalue,
	'minReplicateNum=s'=>\$minReplicateNum,
	'minGroupNum=s'=>\$minGroupNum,
	'outputDir=s'=>\$outputDir,
);

#print join("\t", $taxonId, $jobId, $exptDir, $exptGroupInfo, $foldChange, $pvalue, $minReplicateNum, $minGroupNum, $outputDir) . "\n";
system("rm -rf $outputDir");
system("mkdir -p $outputDir");
# 把exptId读入到group对应的hash中
# root->SRX1203042,SRX1203043,SRX1203045
# leaf->SRX1203049,SRX1203050,SRX1203055
my (@allExptId, %exptToGroup);
my (%groupToExptIdList, %tmpHash, @fieldName, @fieldValue, $line);
open FF, "<$exptGroupInfo";
my $line = <FF>;
chomp($line);
@fieldName = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@fieldValue = split(/\t/, $line);
	for(my $i=0; $i<=$#fieldValue; $i++){
		$tmpHash{$fieldName[$i]} = $fieldValue[$i];
	}

	# 将exptId登记到group中
	$groupToExptIdList{$tmpHash{"groupName"}} .= $tmpHash{"exptId"} . ",";
	# 将exptId登记到exptId总数组中
	$allExptId[$#allExptId+1] = $tmpHash{"exptId"};
	# 建立exptId到group的映射关系
	$exptToGroup{$tmpHash{"exptId"}} = $tmpHash{"groupName"};
}
close FF;

# 获得所有group存入数组，同时去掉exptId列表中末尾的逗号
my (@group, $group);
@group = keys(%groupToExptIdList);
for(my $i=0; $i<=$#group; $i++){
	$groupToExptIdList{$group[$i]} = substr($groupToExptIdList{$group[$i]}, 0, length($groupToExptIdList{$group[$i]}) - 1);
}


# 读取基因的在每个group中tpm，将其按照group存储
my ($geneExprListInGroupHref, %geneExprListInGroup);
$geneExprListInGroupHref=\%geneExprListInGroup;
my $groupName;
# 将所有experiment中gene的readCount读入hash
# geneId  readCount       Coverage        FPKM    TPM
# Os03g0644000    265     7.318612        32.061001       53.069466
my (%geneExptReadCount, $geneExptReadCountHref, $geneId, $exptId, @field, %geneId);
$geneExptReadCountHref = \%geneExptReadCount;
foreach $exptId(@allExptId){

	$groupName = $exptToGroup{$exptId};

	open FF, "<$exptDir/$exptId";
	<FF>;
	while($line=<FF>){
		chomp($line);
		@field = split(/\t/, $line);
		$geneId = $field[0];
		# 将readCount按照geneId -> exptId 映射关系登记下来
		$geneExptReadCountHref->{$geneId}->{$exptId} = $field[1];
		
		# 将geneId登记到总目录中
		$geneId{$geneId} = 1;
		
		# 将tpm按照geneId -> group 映射关系登记下来
		$geneExprListInGroupHref->{$geneId}->{$groupName} .= $field[4] . ",";

	}
	close FF;
}


# 把compare参数分配到threadNum个数组中，以便于生成threadNum个命令
my (@cmdControlPara, @cmdTreatmentPara, @cmdFinalOutputPara, $cmdId);

$cmdId = 0;
my $compareNum=0;
for(my $i=0; $i<=$#group; $i++){
	for(my $j=$i+1; $j<=$#group; $j++){
		print join("\t", $group[$i], $group[$j]) . "\n";
		# ++ 为第cmdId命令准备三个重要参数 ++
		# = 1 controlExptIdList
		$cmdControlPara[$cmdId] .= $groupToExptIdList{$group[$i]} . "_";
		# = 2 treatmentExptIdList
		$cmdTreatmentPara[$cmdId] .= $groupToExptIdList{$group[$j]} . "_";
		# = 3 finalOutputList
		$cmdFinalOutputPara[$cmdId] .= $jobId . "." . $taxonId . "." . $group[$i] . "." . $group[$j] . ".diff.expr.gene.by.DESeq2.tsv" . "_";
		# 命令行编号改变
		$cmdId++;
		$cmdId=0 if($cmdId == $threadNum);
		# 统计compare的数量
		$compareNum++;
	}
}

print "total compare num: $compareNum\n";
print "begin compare!\n";

# 生成启动多个并行运行perl程序的shell脚本文件
my $cmd;
my $launchShellFile = $outputDir . "/$jobId.$taxonId.launch.perlProgam.with.Deseq2.cmd.sh";
open WW, ">$launchShellFile";
for(my $i=0; $i<$threadNum; $i++){

	# 去掉最后的"_"	
	$cmdControlPara[$i] = substr($cmdControlPara[$i], 0, length($cmdControlPara[$i]) - 1);
	$cmdTreatmentPara[$i] = substr($cmdTreatmentPara[$i], 0 , length($cmdTreatmentPara[$i]) - 1);
	$cmdFinalOutputPara[$i]= substr($cmdFinalOutputPara[$i], 0 , length($cmdFinalOutputPara[$i]) - 1);

	$cmd = "nohup perl " . $multiCmpProgram . " " .
		" --taxonId " . $taxonId . 
		" --jobId " . $jobId . 
		" --countDir " . $exptDir . 
		" --controlExptIdListList " . $cmdControlPara[$i] . 
		" --treatmentExptIdListList " . $cmdTreatmentPara[$i] .
		" --diffFoldchange " . $foldChange .
		" --diffQvalue " . $pvalue . 
		" --outputDir " . $outputDir . 
		" --finalOutputList " . $cmdFinalOutputPara[$i] . " 2>/dev/null &\n";
	print WW $cmd;
}
close WW;

# 执行shell脚本，让多个比较并列执行
system("chmod +x $launchShellFile");
system("$launchShellFile");

# 检查比较完成的数量
my $finishedCompareNum;
while(1){
	# 检测完成比较的数量
	$finishedCompareNum = `ls -1 $outputDir/$jobId.$taxonId.*.finished 2>/dev/null |wc -l`;
	chomp($finishedCompareNum);

	print "currentlly finished comparison num: $finishedCompareNum\n";

	if($finishedCompareNum == $compareNum){
		last;
	}else{
		sleep(10);
	}
}

print "finished all deseq2!\n";

my @geneId = keys(%geneId);
my ($groupX, $groupY);
my (%geneGroupCompare, $geneGroupCompareHref);
$geneGroupCompareHref = \%geneGroupCompare;

# 设置循环，列出所有的两两group比较组合
for(my $i=0; $i<=$#group; $i++){

  for(my $j=$i+1; $j<=$#group; $j++){
	
	# 执行DESeq2
	print join("\t", $group[$i], $group[$j]) . "\n";
	my $finalOutput = "$outputDir/$jobId.$taxonId.$group[$i].$group[$j].diff.expr.gene.by.DESeq2.tsv";

	# 结合比较信息提取差异表达的gene存入差异表达统一hash $geneGroupCompareHref
	my @field;
	open FF, "<$finalOutput";
	# [geneId]        "baseMean"        "log2FoldChange"  "lfcSE"            "stat"            "pvalue" "padj"
	# "Os04g0398600"  105651.848306657  8.49836815040118  0.170585767533834  49.8187408789282  0         0
	<FF>;
	while($line = <FF>){
		chomp($line);
		@field = ();
		@field = split(/\t/, $line);
		$geneId = substr($field[0], 1, length($field[0])-2 );
		$groupX = $group[$j];
		$groupY = $group[$i];

		$geneGroupCompareHref->{$geneId}->{$groupX . "/" . $groupY}->{"foldChange"} = 2**($field[2]);
		$geneGroupCompareHref->{$geneId}->{$groupY . "/" . $groupX}->{"foldChange"} = 1/(2**($field[2]));
		$geneGroupCompareHref->{$geneId}->{$groupX . "/" . $groupY}->{"padj"} = $field[6];
		$geneGroupCompareHref->{$geneId}->{$groupY . "/" . $groupX}->{"padj"} = $field[6];
		#<STDIN>;
	}
	close FF;
  
  } # for(my $j
} # for(my $i



open WW, ">$outputDir/$jobId.$taxonId.final.specific.expr.gene.list.by.DESeq2.tsv";
# 对所有基因执行如下操作：
#     逐一检查每个group是否为特异
my ($minPvalue, $maxPvalue, $minFC, $maxFC, $spDirection, $giveUp, $tpmListInAllGroups, $tmpMaxFC, $tmpMinFC, $meanTpm, $exitsCompare);
print WW join("\t", "geneId", "spGroup", "spDirection", "meanTpm", "minPvalue", "maxPvalue", "minFC", "maxFC", "tpmListInAllGroup") . "\n";

foreach $geneId(@geneId){

  # 检查每个groupX是否特异
  foreach $groupX(@group){
    # groupX和其它所有的group比较，是否log2FoldChange都是大于0
    ($minPvalue, $maxPvalue, $minFC, $maxFC, $spDirection, $exitsCompare) = (10, -10, 10000000, -1000000, "", "yes");

    foreach $groupY(@group){

        # 同一group不相比
	next if($groupX eq $groupY);
	
	# 检测是否和groupY存在显著比较结果，如果不存在显著比较结果，那么就放弃
	if(not exists($geneGroupCompareHref->{$geneId}->{$groupX . "/" . $groupY})){
		$exitsCompare = "no";
		last;
	}

	if($minFC > $geneGroupCompareHref->{$geneId}->{$groupX . "/" . $groupY}->{"foldChange"}){
		$minFC = $geneGroupCompareHref->{$geneId}->{$groupX . "/" . $groupY}->{"foldChange"};
	}
	if($maxFC < $geneGroupCompareHref->{$geneId}->{$groupX . "/" . $groupY}->{"foldChange"}){
		$maxFC = $geneGroupCompareHref->{$geneId}->{$groupX . "/" . $groupY}->{"foldChange"};
	}

	if($minPvalue > $geneGroupCompareHref->{$geneId}->{$groupX . "/" . $groupY}->{"padj"}){
		$minPvalue = $geneGroupCompareHref->{$geneId}->{$groupX . "/" . $groupY}->{"padj"};
	}
	if($maxPvalue < $geneGroupCompareHref->{$geneId}->{$groupX . "/" . $groupY}->{"padj"}){
		$maxPvalue = $geneGroupCompareHref->{$geneId}->{$groupX . "/" . $groupY}->{"padj"};
	}
    }

    # 检测groupX是否其他所有group都有显著差异的比较结果
    # 如果无显著比较结果，那么geneId在groupX上不可能显著
    if($exitsCompare eq "no"){
       next;
    }

    # 检测是否groupX比其他所有group都小[即：$maxFC<1($minFC肯定<1)]或者都大[即：$minFC>1($maxFC肯定>1)]
    # 如果不是这样的，那么放弃当前的groupX
    next if(not(($minFC > 1 and $maxFC > 1) or ($minFC < 1 and $maxFC < 1))); 

    #print join("\t", "minFC:" . $minFC, "maxFC:" . $maxFC, "minPvalue:" . $minPvalue, "maxPvalue:" . $maxPvalue);
    #<STDIN>;
    # groupX特异高
    if($minFC > 1 and $maxFC > 1){
    	$spDirection = "high";
    }elsif($minFC < 1 and $maxFC < 1){
        $spDirection = "low";
    }

    # groupY特异低
    if($maxFC < 1){
	# 转换比值为改变的倍数
	if($minFC !=0 ){
		$tmpMaxFC = 1/$minFC;
	}else{
		$tmpMaxFC = 100000000;
	}
	if($maxFC != 0 ){
		$tmpMinFC = 1/$maxFC;
	}else{
		$tmpMinFC = 10000000;
	}
	$maxFC = $tmpMaxFC;
	$minFC = $tmpMinFC;
    }

    # 求出当前groupX中tpm的均值
    my @tpm = ();
    
    @tpm = split(/,/, substr($geneExprListInGroupHref->{$geneId}->{$groupX}, 0, length($geneExprListInGroupHref->{$geneId}->{$groupX})-1) );
    my $meanTpm = sum(@tpm)/($#tpm+1);

    # 将所有group的tpm都列出来
    my $tpmListInAllGroups = "";
    foreach $groupY(@group){
	$tpmListInAllGroups .= $groupY . ":" . substr($geneExprListInGroupHref->{$geneId}->{$groupY}, 0, length($geneExprListInGroupHref->{$geneId}->{$groupY}) - 1) . ";";
    }

    # gene在groupX中特异表达，输出该geneId在groupX及其所有group中表达情况
    print WW join("\t", $geneId, $groupX, $spDirection, sprintf("%.4f", $meanTpm), $minPvalue, $maxPvalue, $minFC, $maxFC, substr($tpmListInAllGroups, 0, length($tpmListInAllGroups))) . "\n";

  } # foreach groupX

}# foreach $geneId

close WW;
