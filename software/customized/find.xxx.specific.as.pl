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
                "--minDltPsi \\\n" .
                "--outputSpecificPsiAs \n";
	exit;
}

my ($specifiedField, $experimentFile, $unitVolum, $minReadCovPerUnitVolum, $psiFileList, 
$minReplicateNum, $minFieldTypeNum, $minDltPsi, $outputSpecificPsiAs);

GetOptions(
        'experimentFile=s'=>\$experimentFile,
	'specifiedField=s'=>\$specifiedField,
        'unitVolum=s'=>\$unitVolum,
        'minReadCovPerUnitVolum=s'=>\$minReadCovPerUnitVolum,
        'psiFileList=s'=>\$psiFileList,
        'minReplicateNum=s'=>\$minReplicateNum,
        'minFieldTypeNum=s'=>\$minFieldTypeNum,
	'minDltPsi=s'=>\$minDltPsi,
	'outputSpecificPsiAs=s'=>\$outputSpecificPsiAs,
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
		$psi=$inclusionNorm/($inclusionNorm+$exclusionNorm);

		# 将Psi值登记到AS的相关fieldType中
		if(not(exists($asToFieldTypeToPsiHref->{$asId}->{$fieldType}))){
			$asToFieldTypeToPsiHref->{$asId}->{$fieldType} = $psi;
		}else{
			$asToFieldTypeToPsiHref->{$asId}->{$fieldType} .= "," . $psi;
		}
	}
	close FF;
}

# 输出包含组织特异的as
my (@asId, @fieldType, $checkReplicate, @psi, $validFieldTypeNum, 
%validFieldTypeToAvgPsi, $validFieldTypeToAvgPsiHref, 
$maxFieldType, $minFieldType, $maxPsi, $minPsi, 
$checkSpecificMax, $checkSpecificMin);
open WW, ">$outputSpecificPsiAs";
print WW join("\t", "ASID", "Regulate", "Tissue/Treatment=Max/MinPsi", "AllFieldTypePsiList") . "\n";
@asId = keys(%asToFieldTypeToPsi);
# 对于某个AS而言：
# (1) 可以既在A组织中特异高psi，也可以在B组织中特异低psi
# (2) 只找到特异高psi的组织，但是没有找到特异低psi的组织
# (3) 只找到特异低psi的组织，但是没有找到特异高psi的组织
# (4) 机没找到特异高psi的组织，也没找到特异低psi的组织
foreach $asId(@asId){
	@fieldType = ();
	@fieldType = keys(%{$asToFieldTypeToPsiHref->{$asId}});

	# 检查每个fieldType中收集的PSI数量是否达到指定值: minReplicateNum
	$validFieldTypeNum = 0;
	%validFieldTypeToAvgPsi = ();
	$validFieldTypeToAvgPsiHref = \%validFieldTypeToAvgPsi;
	# 扫描每一个组织，判断对其做的实验数量是否达到指定值
	foreach $fieldType(@fieldType){
		@psi = ();
		# 获得该组织上做的所有有效实验的psi值
		@psi = split(",", $asToFieldTypeToPsiHref->{$asId}->{$fieldType});
		if($#psi+1 >= $minReplicateNum){
			# 如果该组织上包含的有效实验数量达标，那么将有效组织数增加一，并将psi平均值添加到hash中
			# 否则放弃该组织类型
			$validFieldTypeNum++;
			# 计算该fieldType内的psi平均值，并将平均值登记到hash中
			$validFieldTypeToAvgPsiHref->{$fieldType} = sum(@psi)/($#psi+1);
		}
	}

	# 检查有效fieldType数量是否达到指定值
	next if($validFieldTypeNum < $minFieldTypeNum);

	# 到目前为止，AS达到以下标准：
	# (1) AS在指定数量以上的组织中被检测到；
	# (2) 上述组织中包含的实验重复数也达标
	
	# 接下来，就是找到最高平均PSI值的组织（即特异高PSI组织）和最低平均PSI值的组织
	
	# 
	# 首先找出最大psi平均值和最小psi平均值对应的组织
	$maxPsi = -1;
	$minPsi = 2;
	@fieldType = ();
	@fieldType = keys(%validFieldTypeToAvgPsi);
	foreach $fieldType(@fieldType){
		if($maxPsi < $validFieldTypeToAvgPsiHref->{$fieldType}){
			$maxPsi = $validFieldTypeToAvgPsiHref->{$fieldType};
			$maxFieldType = $fieldType;
		}
		if($minPsi > $validFieldTypeToAvgPsiHref->{$fieldType}){
			$minPsi = $validFieldTypeToAvgPsiHref->{$fieldType};
			$minFieldType = $fieldType;
		}
	}

	# 到目前为止找到了平均最大PSI值的组织maxFieldType和平均最小PSI值的组织minFieldType

	# 检查最大psi平均值是否比其他所有组织的psi值都大0.25
	# 检查最小psi平均值是否比其他所有组织的psi值都小0.25
	$checkSpecificMax = 1;
	$checkSpecificMin = 1;
	foreach $fieldType(@fieldType){
		if($maxFieldType ne $fieldType and $validFieldTypeToAvgPsiHref->{$maxFieldType} - $validFieldTypeToAvgPsiHref->{$fieldType} < $minDltPsi){
			$checkSpecificMax = 0;
		}
		if($minFieldType ne $fieldType and $validFieldTypeToAvgPsiHref->{$fieldType} - $validFieldTypeToAvgPsiHref->{$minFieldType} < $minDltPsi){
			$checkSpecificMin = 0;
		}
	}

	# 如果maxFieldType和minFieldType检查都没通过，那么放弃
	if($checkSpecificMax == 1){
		# 最大psi值对应的组织和最小psi值对应的组织都通过
		# (1) 输出特异高psi的组织及其它组织对应的psi值
		print WW join("\t", $asId, "High", $maxFieldType . "=" . sprintf("%.3f", $validFieldTypeToAvgPsiHref->{$maxFieldType}));
		foreach $fieldType(@fieldType){
			if($fieldType ne $maxFieldType){
				print WW "\t" . $fieldType . "=" . sprintf("%.3f", $validFieldTypeToAvgPsiHref->{$fieldType});
			}
		}
		print WW "\n";
	}
	if($checkSpecificMin == 1){
		print WW join("\t", $asId, "Low", $minFieldType . "=" . sprintf("%.3f", $validFieldTypeToAvgPsiHref->{$minFieldType}));
		foreach $fieldType(@fieldType){
			if($fieldType ne $minFieldType){
				print WW "\t" . $fieldType . "=" . sprintf("%.3f", $validFieldTypeToAvgPsiHref->{$fieldType});
			}
		}
		print WW "\n";

	}
}
close WW;
