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
		"--minReadCoverage 10 \\\n".
                "--changedPsi 0.05 \\\n" .
                "--pvalue 0.05 \\\n" .
                "--minReplicateNum 3\\\n" .
                "--minGroupNum 5 \\\n" .
		"--outputDir specificityByDeseq2Dir\n";
	exit;
}

my ($taxonId, $jobId, $exptDir, $exptGroupInfo, $threadNum, $multiCmpProgram,
$statTest, $foldChange, $pvalue, $changedPsi, $minReadCoverage,
$minReplicateNum, $minGroupNum, $outputDir);

GetOptions(
	'threadNum=s'=>\$threadNum,
	'multiCmpProgram=s'=>\$multiCmpProgram,
        'taxonId=s'=>\$taxonId,
        'jobId=s'=>\$jobId,
        'exptDir=s'=>\$exptDir,
        'exptGroupInfo=s'=>\$exptGroupInfo,
	'statTest=s'=>\$statTest,
	'minReadCoverage=s'=>\$minReadCoverage,
	'changedPsi=s'=>\$changedPsi,
	'pvalue=s'=>\$pvalue,
	'minReplicateNum=s'=>\$minReplicateNum,
	'minGroupNum=s'=>\$minGroupNum,
	'outputDir=s'=>\$outputDir,
);

#print join("\t", $taxonId, $jobId, $exptDir, $exptGroupInfo, $foldChange, $pvalue, $minReplicateNum, $minGroupNum, $outputDir) . "\n";
system("rm -rf $outputDir");
system("mkdir -p $outputDir");
# (1) groupToExptIdList
# root->SRX1203042,SRX1203043,SRX1203045
# (2) exptToGroup
# SRX1203042 -> root
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

	# 将exptIdList登记到group中
	$groupToExptIdList{$tmpHash{"groupName"}} .= $tmpHash{"exptId"} . ",";

	# 将exptId登记到exptId总数组中
	$allExptId[$#allExptId+1] = $tmpHash{"exptId"};

	# 建立exptId到group的映射关系
	$exptToGroup{$tmpHash{"exptId"}} = $tmpHash{"groupName"};
}
close FF;

# 获得所有group
my (@group, $group);
@group = keys(%groupToExptIdList);

# 处理group指向exptIdList，去掉末尾的","
for(my $i=0; $i<=$#group; $i++){
	$groupToExptIdList{$group[$i]} = substr($groupToExptIdList{$group[$i]}, 0, length($groupToExptIdList{$group[$i]}) - 1);
}

# 读取AS在每个group中psi列表，将其按照group存储
my ($asPsiListInGroupHref, %asPsiListInGroup);
$asPsiListInGroupHref=\%asPsiListInGroup;
my $groupName;
# ASID    		IJC_SAMPLE_1    SJC_SAMPLE_1    IncFormLen      SkipFormLen     Experiment
# OSJGA3SS0000000001    0       	6       	103     	99      	DRX103903
my (%asExptPsi, $asExptPsiHref, $asId, $exptId, @field, %asId, $IJC_SAMPLE, $SJC_SAMPLE, $IncFormLen, $SkipFormLen, $expermentId, $psi);
$asExptPsiHref = \%asExptPsi;
foreach $exptId(@allExptId){

	# 获得当前expt对应的groupName	
	$groupName = $exptToGroup{$exptId};

	open FF, "<$exptDir/$exptId";
	<FF>;
	while($line=<FF>){

		chomp($line);

		($asId, $IJC_SAMPLE, $SJC_SAMPLE, $IncFormLen, $SkipFormLen, $expermentId) = split(/\t/, $line);

		# 检测read覆盖数是否达标
		next if($IJC_SAMPLE + $SJC_SAMPLE < $minReadCoverage);

		# 计算psi值
		$psi = ($IJC_SAMPLE/$IncFormLen)/($IJC_SAMPLE/$IncFormLen + $SJC_SAMPLE/$SkipFormLen);

		# 将psi按照asId -> exptId 映射关系登记下来
		$asExptPsiHref->{$asId}->{$exptId} = $psi;
		
		# 将asId登记到总目录中
		$asId{$asId} = 1;
		
		# 将psi登记到相应的group中
		$asPsiListInGroupHref->{$asId}->{$groupName} .= $psi . ",";
	}
	close FF;
}

#  asPsiListInGroupHref: 用于求每个group中psi值的平均值
#  @group :              用于登记所有group
#  groupToExptIdList:    用于列出比较时对应的exptId列表

# 双重循环列出所有比较
# 把compare参数分配到threadNum个数组中，以便于生成threadNum个命令
my (@cmdControlPara, @cmdTreatmentPara, @cmdFinalOutputPara, $cmdId);
$cmdId = 0;
my $compareNum=0;
for(my $i=0; $i<=$#group; $i++){
	for(my $j=$i+1; $j<=$#group; $j++){
		# print join("\t", $group[$i], $group[$j]) . "\n";
		# ++ 为第cmdId命令准备三个重要参数 ++
		# = 1 controlExptIdList
		$cmdControlPara[$cmdId] .= $groupToExptIdList{$group[$i]} . "_";
		# = 2 treatmentExptIdList
		$cmdTreatmentPara[$cmdId] .= $groupToExptIdList{$group[$j]} . "_";
		# = 3 finalOutputList
		$cmdFinalOutputPara[$cmdId] .= $jobId . "." . $taxonId . "." . $group[$i] . "." . $group[$j] . ".diff.splicing.as.by.psi.tsv" . "_";
		# 命令行编号改变
		$cmdId++;
		$cmdId=0 if($cmdId == $threadNum);
		# 统计compare的数量
		$compareNum++;
	}
}

# print "total compare num: $compareNum\n";

# 生成启动多个并行运行perl程序的shell脚本文件
my $cmd;
my $launchShellFile = $outputDir . "/$jobId.$taxonId.launch.perlProgam.with.psi.cmd.sh";
open WW, ">$launchShellFile";
for(my $i=0; $i<$threadNum; $i++){

	# 去掉最后的"_"	
	$cmdControlPara[$i] = substr($cmdControlPara[$i], 0, length($cmdControlPara[$i]) - 1);
	$cmdTreatmentPara[$i] = substr($cmdTreatmentPara[$i], 0 , length($cmdTreatmentPara[$i]) - 1);
	$cmdFinalOutputPara[$i]= substr($cmdFinalOutputPara[$i], 0 , length($cmdFinalOutputPara[$i]) - 1);

	$cmd = "nohup perl " . $multiCmpProgram . " " .
		" --taxonId " . $taxonId . 
		" --jobId " . $jobId . 
		" --psiDir " . $exptDir . 
		" --controlExptIdListList " . $cmdControlPara[$i] . 
		" --treatmentExptIdListList " . $cmdTreatmentPara[$i] .
		" --changedPsi " . $changedPsi .
		" --minReadCoverage " . $minReadCoverage .
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

	# print "currentlly finished comparison num: $finishedCompareNum\n";

	if($finishedCompareNum == $compareNum){
		last;
	}else{
		sleep(10);
	}
}

# print "finished!\n";

my @asId = keys(%asId);
my ($groupX, $groupY);
my (%asGroupCompare, $asGroupCompareHref);
$asGroupCompareHref = \%asGroupCompare;

# 设置循环，列出所有的两两group比较组合
for(my $i=0; $i<=$#group; $i++){
  for(my $j=$i+1; $j<=$#group; $j++){
	
	#print join("\t", $group[$i], $group[$j]) . "\n";
	my $finalOutput = "$outputDir/$jobId.$taxonId.$group[$i].$group[$j].diff.splicing.as.by.psi.tsv";

	# 结合比较信息提取差异AS存入hash
	my @field;
	open FF, "<$finalOutput";
	# as_Id            control         treatment              changedPsi              pvalue
	# OSJGRI0000008419 0.12,0.10,0.0   76.4,9.04,6.07,90.12   15951.6176988758        0.00290610522289939
	<FF>;
	while($line = <FF>){
		chomp($line);
		@field = ();
		@field = split(/\t/, $line);
		$asId = $field[0];

		$groupX = $group[$j];
		$groupY = $group[$i];

		$asGroupCompareHref->{$asId}->{$groupX . "/" . $groupY}->{"changed"} = $field[3];
		$asGroupCompareHref->{$asId}->{$groupY . "/" . $groupX}->{"changed"} = -1 * $field[3];

		$asGroupCompareHref->{$asId}->{$groupX . "/" . $groupY}->{"padj"} = $field[4];
		$asGroupCompareHref->{$asId}->{$groupY . "/" . $groupX}->{"padj"} = $field[4];
		#<STDIN>;
	}
	close FF;
  
  } # for(my $j
} # for(my $i


open WW, ">$outputDir/$jobId.$taxonId.final.specific.splicing.as.list.by.psi.tsv";

# 逐一检查每个group是否为特异
my ($minPvalue, $maxPvalue, $minCP, $maxCP, $spDirection, $giveUp, $psiListInAllGroups, $tmpMaxCP, $tmpMinCF, $meanPsi, $exitsCompare);
print WW join("\t", "asId", "spGroup", "spDirection", "meanPsi", "minPvalue", "maxPvalue", "minCP", "maxCP", "psiListInAllGroup") . "\n";

foreach $asId(@asId){

  # 检查每个groupX是否特异
  foreach $groupX(@group){
    # groupX和其它所有的group比较，是否log2FoldChange都是大于0
    ($minPvalue, $maxPvalue, $minCP, $maxCP, $spDirection, $exitsCompare) = (10, -10, 10000000, -1000000, "", "yes");
    # 扫描所有其他groupY
    foreach $groupY(@group){
        # 同一group不相比
	next if($groupX eq $groupY);

	# 如果asPsiListInGroupHref->{$asId}->{$groupY}不存在，则表示在groupY上未检测到AS
	# 此时，如果$asGroupCompareHref->{$asId}->{$groupX . "/" . $groupY}不存在，那么并不表示在groupX对groupY不显著
	# 因此：只有$asPsiListInGroupHref->{$asId}->{$groupY}存在，且$groupX . "/" . $groupY不存在时才放弃
	#             没有检测到显著比较结果                                              groupY上存在PSI
	if(not exists($asGroupCompareHref->{$asId}->{$groupX . "/" . $groupY}) and exists($asPsiListInGroupHref->{$asId}->{$groupY})){
		$exitsCompare = "no";
		last;
	}

	if($minCP > $asGroupCompareHref->{$asId}->{$groupX . "/" . $groupY}->{"changed"}){
		$minCP = $asGroupCompareHref->{$asId}->{$groupX . "/" . $groupY}->{"changed"};
	}

	if($maxCP < $asGroupCompareHref->{$asId}->{$groupX . "/" . $groupY}->{"changed"}){
		$maxCP = $asGroupCompareHref->{$asId}->{$groupX . "/" . $groupY}->{"changed"};
	}

	if($minPvalue > $asGroupCompareHref->{$asId}->{$groupX . "/" . $groupY}->{"padj"}){
		$minPvalue = $asGroupCompareHref->{$asId}->{$groupX . "/" . $groupY}->{"padj"};
	}
	if($maxPvalue < $asGroupCompareHref->{$asId}->{$groupX . "/" . $groupY}->{"padj"}){
		$maxPvalue = $asGroupCompareHref->{$asId}->{$groupX . "/" . $groupY}->{"padj"};
	}
    }

    # 检测groupX是否其他所有group都有显著差异的比较结果. 如果无显著比较结果，那么asId在groupX上不可能显著
    if($exitsCompare eq "no"){
       next;
    }

    # 检测是否groupX比其他所有group都小[即：$maxFC<1($minFC肯定<1)]或者都大[即：$minFC>1($maxFC肯定>1)]
    # 如果不是这样的，那么放弃当前的groupX
    next if(not(($minCP > 0 and $maxCP > 0) or ($minCP < 0 and $maxCP < 0))); 

    #print join("\t", "minFC:" . $minFC, "maxFC:" . $maxFC, "minPvalue:" . $minPvalue, "maxPvalue:" . $maxPvalue);
    #<STDIN>;
    # groupX特异高
    if($minCP > 0 and $maxCP > 0){
    	$spDirection = "high";
    }elsif($minCP < 0 and $maxCP < 0){
        $spDirection = "low";
    }

    # groupY特异低
    if($maxCP < 0){
	($minCP, $maxCP) = (-1*$maxCP, -1*$minCP);
    }

    # 求出当前groupX中psi的均值
    my @psi = ();
    
    @psi = split(/,/, substr($asPsiListInGroupHref->{$asId}->{$groupX}, 0, length($asPsiListInGroupHref->{$asId}->{$groupX})-1) );
    my $meanTpm = sum(@psi)/($#psi+1);

    # 将所有group的psi都列出来
    my $psiListInAllGroups = "";
    foreach $groupY(@group){
	if(exists($asPsiListInGroupHref->{$asId}->{$groupY})){
		$psiListInAllGroups .= $groupY . ":" . substr($asPsiListInGroupHref->{$asId}->{$groupY}, 0, length($asPsiListInGroupHref->{$asId}->{$groupY}) - 1) . ";";
	}
    }

    # gene在groupX中特异表达，输出该asId在groupX及其所有group中表达情况
    print WW join("\t", $asId, $groupX, $spDirection, sprintf("%.4f", $meanPsi), $minPvalue, $maxPvalue, $minCP, $maxCP, substr($psiListInAllGroups, 0, length($psiListInAllGroups))) . "\n";

  } # foreach groupX

}# foreach $asId

close WW;
