#!/usr/bin/perl
use strict;
use List::Util qw/max min sum maxstr minstr shuffle/;
use Getopt::Long;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--specifiedField TreatmentType \\\n" .
                "--experimentFile \\\n" .
		"--unitVolum \\\n" .
		"--minReadCovPerUnitVolum \\\n" .
                "--psiFileList \\\n" .
                "--minReplicateNum \\\n" .
                "--minFieldTypeNum \\\n" .
                "--outputCovMatrix \\\n" .
		"--outputConstantAs \n";
	exit;
}

my ($specifiedField, $experimentFile, $unitVolum, $minReadCovPerUnitVolum, $psiFileList, 
$minReplicateNum, $minFieldTypeNum, $minDltPsi, $outputCovMatrix, $outputConstantAs);

GetOptions(
        'experimentFile=s'=>\$experimentFile,
	'specifiedField=s'=>\$specifiedField,
        'unitVolum=s'=>\$unitVolum,
        'minReadCovPerUnitVolum=s'=>\$minReadCovPerUnitVolum,
        'psiFileList=s'=>\$psiFileList,
        'minReplicateNum=s'=>\$minReplicateNum,
        'minFieldTypeNum=s'=>\$minFieldTypeNum,
	'outputCovMatrix=s'=>\$outputCovMatrix,
	'outputConstantAs=s'=>\$outputConstantAs,
);
my ($line, $i, @nameField, @valueField);
##  == 1 获得experiment的组织信息 ===
# 将experiment读入hash，可以通过experimentId直接获得指定field类型
my (%exprimentInfo, $experimentInfoHref, %tmpHash, $exptId);
$experimentInfoHref = \%exprimentInfo;
open FF, "<$experimentFile";
$line = <FF>;
chomp($line);
@nameField = split(/\t/, $line);
while($line=<FF>){
	chomp($line);
	@valueField = split(/\t/, $line);
	for($i=0; $i<=$#valueField; $i++){
		$tmpHash{$nameField[$i]} = $valueField[$i];
	}
	$experimentInfoHref->{$tmpHash{"Experiment"}}->{$specifiedField} = $tmpHash{$specifiedField};
	$experimentInfoHref->{$tmpHash{"Experiment"}}->{"minReadCoverage"} = $tmpHash{"mappedBases"}/$unitVolum * $minReadCovPerUnitVolum;
}
close FF;


## == 2 构建as->fieldType->psiList hash ===
# ASID    IJC_SAMPLE_1    SJC_SAMPLE_1    IncFormLen      SkipFormLen     Experiment
# ATHAA3SS0000000001      549     189     398     100     SRX485073
my (%asToFieldTypeToPsi, $asToFieldTypeToPsiHref, $asId, $experimentId, $fieldType);
my ($inclusionNorm, $exclusionNorm, $psi);
my (@psiFile, $psiFile);
$asToFieldTypeToPsiHref = \%asToFieldTypeToPsi;
@nameField = ();
@valueField = ();
%tmpHash = ();
@psiFile = split(/,/, $psiFileList);
foreach $psiFile(@psiFile){
	open FF, "<$psiFile";
	$line = <FF>;
	chomp($line);
	@nameField = split(/\t/, $line);
	while($line=<FF>){
		chomp($line);
		@valueField = split(/\t/, $line);
		for($i=0; $i<=$#valueField; $i++){
			$tmpHash{$nameField[$i]} = $valueField[$i];
		}

		# 如果该experiment不被关注，那么放弃
		next if(not(exists($experimentInfoHref->{$tmpHash{"Experiment"}})));
		# 如果inclusion和exclusion都没有检测到，那么放弃
		next if($tmpHash{"IJC_SAMPLE_1"} == 0 and $tmpHash{"SJC_SAMPLE_1"} == 0);

		# 如果inclusion和exclusion都没有达到最低覆盖度，那么表示该AS没有被有效检测，也放弃
		next if($tmpHash{"IJC_SAMPLE_1"} + $tmpHash{"SJC_SAMPLE_1"} < $experimentInfoHref->{$tmpHash{"Experiment"}}->{"minReadCoverage"});

		# 计算PSI值
		$experimentId = $tmpHash{"Experiment"};
		$fieldType = $experimentInfoHref->{$experimentId}->{$specifiedField};
		$asId = $tmpHash{"ASID"};

		$inclusionNorm = $tmpHash{"IJC_SAMPLE_1"}/$tmpHash{"IncFormLen"};
		$exclusionNorm = $tmpHash{"SJC_SAMPLE_1"}/$tmpHash{"SkipFormLen"};
		$psi=sprintf("%.3f", $inclusionNorm/($inclusionNorm+$exclusionNorm));

		# 将Psi值登记到AS的相关fieldType中
		if(not(exists($asToFieldTypeToPsiHref->{$asId}->{$fieldType}))){
			$asToFieldTypeToPsiHref->{$asId}->{$fieldType} = $psi;
		}else{
			$asToFieldTypeToPsiHref->{$asId}->{$fieldType} .= "," . $psi;
		}
	}
	close FF;
}

# 1、列出符合基本条件的AS:(1)组织类型数量要达标;(2)每个组织内样品数量
# 2、将这些AS的所有psi值列出来，然后计算协方差
# 	SE0001	0.676,0.876,0.345 	0.23
# 3、根据所有AS的协方差值的分布，找到1/4处的阈值
my (@asId, @fieldType, $checkReplicate, @psi, $validFieldTypeNum, 
%validFieldTypeToPsiList, $validFieldTypeToPsiListHref, 
$psiList);

open COV, ">$outputCovMatrix";
@asId = keys(%asToFieldTypeToPsi);
# 对于某个AS而言：
# (1) 可以既在A组织中特异高psi，也可以在B组织中特异低psi
# (2) 只找到特异高psi的组织，但是没有找到特异低psi的组织
# (3) 只找到特异低psi的组织，但是没有找到特异高psi的组织
# (4) 机没找到特异高psi的组织，也没找到特异低psi的组织
foreach $asId(@asId){
	# 获得asId关联的所有组织类型
	@fieldType = ();
	@fieldType = keys(%{$asToFieldTypeToPsiHref->{$asId}});

	# 检查每个fieldType中收集的PSI数量是否达到指定值: minReplicateNum
	$validFieldTypeNum = 0;
	%validFieldTypeToPsiList = ();
	$validFieldTypeToPsiListHref = \%validFieldTypeToPsiList;

	# 扫描每一个组织，判断对其做的实验数量是否达到指定值
	foreach $fieldType(@fieldType){
		@psi = ();
		# 获得该组织上做的所有有效实验的数量
		@psi = split(",", $asToFieldTypeToPsiHref->{$asId}->{$fieldType});
		if($#psi+1 >= $minReplicateNum){
			# 如果该组织上包含的有效实验数量达标，那么将有效组织数增加1，并将该组织的psiList登记到validFieldTypeToPsiList中
			# 否则放弃该组织
			$validFieldTypeNum++;
			$validFieldTypeToPsiListHref->{$fieldType} = $asToFieldTypeToPsiHref->{$asId}->{$fieldType};
		}
	}

	# 检查有效fieldType数量是否达到指定值
	next if($validFieldTypeNum < $minFieldTypeNum);

	# 到目前为止，AS达到以下标准：
	# (1) AS在指定数量以上的组织中被检测到；
	# (2) 上述组织中包含的实验重复数也达标
	
	# 计算当前AS的psi值的COV值	
	@fieldType = ();
	@fieldType = keys(%validFieldTypeToPsiList);
	$psiList = "";
	foreach $fieldType(@fieldType){
		if($psiList eq ""){
			$psiList = $validFieldTypeToPsiListHref->{$fieldType};
		}else{
			$psiList .= "," . $validFieldTypeToPsiListHref->{$fieldType};
		}		
	}

	# 计算变异系数cov
	my ($psiAvg, $cov) = (-1, -1);
	&calCov($psiList, \$psiAvg, \$cov);
	print COV join("\t", $asId, $psiList, $cov, $psiAvg) . "\n";
}
close COV;

my $totalAsNum=`cat $outputCovMatrix |wc -l`;
my $quarterNum=int($totalAsNum/4);
my $cmd = "echo -e \"asId\tCV\tAvgPsi\" > $outputConstantAs";
system($cmd);
$cmd = "sort -k3,3n $outputCovMatrix |awk -F \'\\t\' \'{print \$1\"\\t\"\$3\"\\t\"\$4}\' | head -n $quarterNum >> $outputConstantAs";
system($cmd);

# 计算变异系数
sub calCov{
	my ($psiList,  $avgPsi, $cov) = @_;
	my (@psi, $psi, $psiU, $sum, $psiVar);
	# 求均值
	@psi = split(/,/, $psiList);
	$psiU = (sum(@psi))/($#psi+1);
	$$avgPsi = $psiU;
	if($psiU == 0){
		$$cov = 0;
		return;
	}
	# 求方差
	$sum = 0;
	foreach $psi(@psi){
		$sum+=($psi-$psiU)*($psi-$psiU)
	}
	$psiVar = ($sum/($#psi))**0.5;
	# 返回变异系数
	$$cov = $psiVar/$psiU;
}



