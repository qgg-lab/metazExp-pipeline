#!/usr/bin/perl
use strict;
use Statistics::R;
use List::Util qw/max min sum maxstr minstr shuffle/;
use Getopt::Long;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--specifiedField TreatmentType \\\n" .
		"--ignoreTypeList seedling \\\n" .
		"--minDltPsi 2\\\n" .
                "--experimentFile \\\n" .
		"--unitVolum \\\n" .
		"--minReadCovPerUnitVolum \\\n" .
                "--psiFileList \\\n" .
                "--minReplicateNum \\\n" .
                "--minFieldTypeNum \\\n" .
                "--outputSpecificPsiAs \n";
	exit;
}

my ($specifiedField, $experimentFile, $unitVolum, $minReadCovPerUnitVolum, $psiFileList, $ignoreTypeList, $minDltPsi,
$minReplicateNum, $minFieldTypeNum, $outputSpecificPsiAs);

GetOptions(
        'experimentFile=s'=>\$experimentFile,
	'ignoreTypeList=s'=>\$ignoreTypeList,
	'minDltPsi=s'=>\$minDltPsi,
	'specifiedField=s'=>\$specifiedField,
        'unitVolum=s'=>\$unitVolum,
        'minReadCovPerUnitVolum=s'=>\$minReadCovPerUnitVolum,
        'psiFileList=s'=>\$psiFileList,
        'minReplicateNum=s'=>\$minReplicateNum,
        'minFieldTypeNum=s'=>\$minFieldTypeNum,
	'outputSpecificPsiAs=s'=>\$outputSpecificPsiAs,
);
my ($line, $i, @nameField, @valueField);
##  == 1 获得experiment的组织信息 ===
# 将experiment读入hash，可以通过experimentId直接获得指定field类型
# $experimentInfoHref->{"ERX12345"}->{"Tissue"} = "leaf";
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

	# 如果该组织/处理不是关注的类型，那么放弃登记
	next if(index($ignoreTypeList, $tmpHash{$specifiedField})>=0);

	# 登记
	$experimentInfoHref->{$tmpHash{"Experiment"}}->{$specifiedField} = $tmpHash{$specifiedField};
	$experimentInfoHref->{$tmpHash{"Experiment"}}->{"minReadCoverage"} = $tmpHash{"Base"}/$unitVolum * $minReadCovPerUnitVolum;
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

		# 如果该experiment没有登记组织类型则不被关注，那么放弃
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
		# $asToFieldTypeToPsiHref->{"ATHAA3SS0000002300"}->{"leaf"} = "0.23,0.89,0.45,012";
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
print WW join("\t", "pvalue0.01", "pvalue0.03", "pvalue0.03", "ASID", "Regulate", "$specifiedField=PSI|exptNum|minDltPsi", "1stTissue|avgPSI|exptNum|pvalue", "1stTissue|avgPSI|exptNum|pvalue", "……") . "\n";
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
	# $asToFieldTypeToPsiHref->{"ATHAA3SS0000002300"}->{"leaf"} = "0.23,0.89,0.45,0.12";
	# $asToFieldTypeToPsiHref->{"ATHAA3SS0000002300"}->{"root"} = "0.28,0.93,0.15,0.33";
	
	# 收集样本重复数（即PSI值数量）达到规定数的组织，然后计算PSI值的平均值登记到validFieldTypeToAvgPsi
	$validFieldTypeNum = 0;
	%validFieldTypeToAvgPsi = ();
	$validFieldTypeToAvgPsiHref = \%validFieldTypeToAvgPsi;
	my $validFieldTypeList = ""; #登记样本重复数符合条件的组织数量
	foreach $fieldType(@fieldType){
		@psi = ();
		# 获得该组织上做的所有有效实验的psi值
		@psi = split(",", $asToFieldTypeToPsiHref->{$asId}->{$fieldType});
		if($#psi+1 >= $minReplicateNum){
			# 如果该组织上包含的有效实验数量达标，那么将有效组织数增加一，并将psi平均值添加到hash中
			# 否则放弃该组织类型
			$validFieldTypeNum++;
			$validFieldTypeList .= $fieldType . ",";
			# 计算该fieldType内的psi平均值，并将平均值登记到hash中
			# $validFieldTypeToAvgPsiHref->{"leaf"} = 0.23
			# $validFieldTypeToAvgPsiHref->{"root"} = 0.14
			# ……
			$validFieldTypeToAvgPsiHref->{$fieldType} = sum(@psi)/($#psi+1);
		}
	}
	$validFieldTypeList = substr($validFieldTypeList, 0, length($validFieldTypeList) - 1);


	# 检查有效当前asId下有小组织数量是否达到指定值
	next if($validFieldTypeNum < $minFieldTypeNum);

	# 到目前为止，AS达到以下标准：
	# (1) AS在指定数量以上的组织中被检测到；
	# (2) 上述组织中包含的实验重复数也达标
	

	# 找出最大psi平均值和最小psi平均值对应的组织
	$maxPsi = -1;
	$minPsi = 2;
	@fieldType = ();
	@fieldType = keys(%validFieldTypeToAvgPsi);
	# $asToFieldTypeToPsiHref->{"ATHAA3SS0000002300"}->{"leaf"} = "0.23,0.89,0.45,0.12";
	# $asToFieldTypeToPsiHref->{"ATHAA3SS0000002300"}->{"root"} = "0.28,0.93,0.15,0.33";
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


########## 用t检验检查最大值的组织特异性 ###############
	# 检查最大值Psi和其它psi值之差是否都小于$minDltPsi
	my $highEnoughTag = 1;
	my $highMinDltPsi = 1;
	foreach $fieldType(@fieldType){
		if($maxPsi-$validFieldTypeToAvgPsiHref->{$fieldType}< $minDltPsi and $fieldType ne $maxFieldType){
			$highEnoughTag = 0;
		}
		if($fieldType ne $maxFieldType){
			$highMinDltPsi = $maxPsi-$validFieldTypeToAvgPsiHref->{$fieldType} if($maxPsi-$validFieldTypeToAvgPsiHref->{$fieldType} < $highMinDltPsi);
		}
	}


	if($highEnoughTag == 1){
		my($rltListBy3Cutoff, $pvalueList) = ("N\tN\tN", "");
		# $checkRltListBy3Cutoff = "N\tY\tY"
		# $pvalueList = "0.02\t0.12\t0.15\t0.13\t0.06"
		&tTestSpecific($asToFieldTypeToPsiHref, $validFieldTypeToAvgPsiHref, $asId, $validFieldTypeList, $maxFieldType, "high", \$rltListBy3Cutoff, \$pvalueList);
		# asToFieldTypeToPsiHref->{$asId}->{"..."} 存放的是所有asId下，所有组织下的psi值列表
	        my @tmpArr = ();
	        @tmpArr = split(/,/, $asToFieldTypeToPsiHref->{$asId}->{$maxFieldType});
	        my $exptNum = $#tmpArr + 1;
		print WW join("\t", $rltListBy3Cutoff, $asId, "High", $maxFieldType . "=" . sprintf("%.3f", $validFieldTypeToAvgPsiHref->{$maxFieldType}) . "|" . $exptNum . "|" . sprintf("%.3f", $highMinDltPsi));
		@fieldType = ();
		@fieldType = split(/,/, $validFieldTypeList);
		my $pvalue = "";
		my @pvalue=();
		@pvalue = split(/,/, $pvalueList);
		my $pk = 0;
		#"1stTissue|avgPSI|exptNum|pvalue"
		for(my $m=0; $m<=$#fieldType;  $m++){
			$pvalue = $pvalue[$pk];
			$pk++;
			$fieldType = $fieldType[$m];
			if($fieldType eq $maxFieldType){
				$pk--;
				next;
			}
			my @tt = ();
			@tt = split(/,/, $asToFieldTypeToPsiHref->{$asId}->{$fieldType}); # 用于统计psi值的个数
			my $exptNum = $#tt+1;
			print WW "\t" . join("|", $fieldType, sprintf("%.3f",$validFieldTypeToAvgPsiHref->{$fieldType}), $exptNum, $pvalue);
		}
		print WW "\n";
	}





######### 用t检验检查最小值的组织特异性 #################
	# 检查最小值Psi和其它psi值之差是否都小于$minDltPsi
	my $lowEnoughTag = 1;
	my $lowMinDltPsi = 1;
	foreach $fieldType(@fieldType){
		if($validFieldTypeToAvgPsiHref->{$fieldType} - $minPsi < $minDltPsi and $fieldType ne $minFieldType){
			$lowEnoughTag = 0;
		}
		if($fieldType ne $minFieldType){
			$lowMinDltPsi = $validFieldTypeToAvgPsiHref->{$fieldType} - $minPsi if($validFieldTypeToAvgPsiHref->{$fieldType} - $minPsi<$lowMinDltPsi);
		}
	}

	if($lowEnoughTag == 1){
		my($rltListBy3Cutoff, $pvalueList) = ("N\tN\tN", "");
		# $checkRltListBy3Cutoff = "N\tY\tY"
		# $pvalueList = "0.02\t0.12\t0.15\t0.13\t0.06"
		&tTestSpecific($asToFieldTypeToPsiHref, $validFieldTypeToAvgPsiHref, $asId, $validFieldTypeList, $minFieldType, "low", \$rltListBy3Cutoff, \$pvalueList);
		# asToFieldTypeToPsiHref->{$asId}->{"..."} 存放的是所有asId下，所有组织下的psi值列表
		my @tmpArr = ();
		@tmpArr = split(/,/, $asToFieldTypeToPsiHref->{$asId}->{$minFieldType});
		my $exptNum = $#tmpArr + 1;
		print WW join("\t", $rltListBy3Cutoff, $asId, "Low", $minFieldType . "=" . sprintf("%.3f", $validFieldTypeToAvgPsiHref->{$minFieldType}) . "|" . $exptNum . "|" . sprintf("%.3f", $lowMinDltPsi));
		@fieldType = ();
		@fieldType = split(/,/, $validFieldTypeList);
	
		my $pvalue = "";
		my @pvalue=();
		@pvalue = split(/,/, $pvalueList);
		#"1stTissue|avgPSI|exptNum|pvalue"
		my $pk = 0;
		for(my $m=0; $m<=$#fieldType;  $m++){
			$pvalue = $pvalue[$pk];
			$pk++;
			$fieldType = $fieldType[$m];
			if($fieldType eq $minFieldType){
				$pk--;
				next;
			}
			my @tt = ();
			@tt = split(/,/, $asToFieldTypeToPsiHref->{$asId}->{$fieldType}); # 用于统计psi值的个数
			my $exptNum = $#tt+1;
			print WW "\t" . join("|", $fieldType, sprintf("%.3f",$validFieldTypeToAvgPsiHref->{$fieldType}), $exptNum, $pvalue);
		}
		print WW "\n";
	}
}
close WW;

# 
sub tTestSpecific{
	my ($asToFieldTypeToPsiHref, $avgHref, $asId, $validFieldTypeList, $specificFieldType, $regulateType, $rltListBy3Cutoff, $pvalueList) = @_;
	$$rltListBy3Cutoff = "";
	$$pvalueList = "";
	my ($pvalue);
	my @fieldType = split(/,/, $validFieldTypeList);
	my ($pvalue001Flag, $pvalue003Flag, $pvalue005Flag) = ("Y", "Y", "Y");
	my $specificFieldPsiList = $asToFieldTypeToPsiHref->{$asId}->{$specificFieldType};
	my $R = Statistics::R->new( shared => 1);

	#print join("\t", "validFieldTypeList:" . $validFieldTypeList, "specificFieldType:" . $specificFieldType);
	#<STDIN>;
	foreach my $fieldType(@fieldType){

		next if($fieldType eq $specificFieldType);

		my $psiList = $asToFieldTypeToPsiHref->{$asId}->{$fieldType};

my $flag = "greater";
if(uc($regulateType) eq "LOW"){
	$flag = "less";
}
my $cmd=<<EOF;
x<-c($specificFieldPsiList)
y<-c($psiList)
pvalue<-t.test(x, y, alternative = "$flag", paired=FALSE)\$p.value
EOF
		#print $cmd;
		#<STDIN>;
		# 两组数据内部相同，且组间相同，则表示无差异，强制设置pvalue = 1
		if(&innerSame($specificFieldPsiList)==1 and &innerSame($psiList)==1 and $avgHref->{$specificFieldType} == $avgHref->{$fieldType}){ 
			$pvalue = 1;
		# 两组数据内部相同，但是组间不相同，则表示两个不相同的常数，强制设置pvalue=0
		}elsif(&innerSame($specificFieldPsiList)==1 and &innerSame($psiList)==1 and $avgHref->{$specificFieldType} != $avgHref->{$fieldType}){
			$pvalue = 0;
		# 两组数据只要有1组内部不相同，那么就可以执行t检验
		}elsif(&innerSame($specificFieldPsiList)!=1 or &innerSame($psiList)!=1){
			$R->run($cmd);
			$pvalue=$R->get('pvalue');
		}else{
			$R->run($cmd);
			$pvalue=$R->get('pvalue');
		}
		$pvalue001Flag = "N" if($pvalue > 0.01);
		$pvalue003Flag = "N" if($pvalue > 0.03);
		$pvalue005Flag = "N" if($pvalue > 0.05);

		#  检测倍数是否符合要求
		
		if($$pvalueList eq ""){
			$$pvalueList = sprintf("%.3f", $pvalue);
		}else{
			$$pvalueList .= "," . sprintf("%.3f", $pvalue);
		}
	}
	$$rltListBy3Cutoff = join("\t", $pvalue001Flag, $pvalue003Flag, $pvalue005Flag);
	#print join("\t", "rltListBy3Cutoff:" . $$rltListBy3Cutoff, "pvalueList:" . $$pvalueList);
	#<STDIN>;
}



# check innerPsi
sub innerSame{
        my ($psiList) = @_;
        my (%psi, $psi, @psi);
        @psi = split(/,/, $psiList);
        foreach $psi(@psi){
                $psi{$psi} = 1;
        }
        @psi = ();
        @psi = keys(%psi);
        if($#psi == 0){
                return 1;
        }else{
                return 0;
        }
}
