#!/usr/bin/perl
use strict;
use Statistics::R;
use List::Util qw/max min sum maxstr minstr shuffle/;
use Getopt::Long;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--geneExpressTsv  \\\n" .
                "--experimentFile \\\n" .
		"--specifiedField Tissue\\\n" .
                "--minReplicateNum \\\n" .
                "--minFieldTypeNum \\\n" .
		"--outputExprCVMatrix \\\n" .
		"--outputConstant \n";
	exit;
}

my ($specifiedField, $experimentFile, $geneExpressTsv,
$minReplicateNum, $minFieldTypeNum,
$outputExprCVMatrix, $outputConstant);

GetOptions(
        'experimentFile=s'=>\$experimentFile,
	'specifiedField=s'=>\$specifiedField,
        'geneExpressTsv=s'=>\$geneExpressTsv,
        'minReplicateNum=s'=>\$minReplicateNum,
        'minFieldTypeNum=s'=>\$minFieldTypeNum,
	'outputExprCVMatrix=s'=>\$outputExprCVMatrix,
	'outputConstant=s'=>\$outputConstant,
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
#	print join("\t", "ExperimentId:" . $tmpHash{"Experiment"}, $specifiedField . ":" . $tmpHash{$specifiedField}) . "\n";
#	print "experiment->tissue:" . $experimentInfoHref->{$tmpHash{"Experiment"}}->{$specifiedField} . "\n";
#	<STDIN>;
}
close FF;


## == 2 构建gene->fieldType->tpmList ===
# GeneId, ExptId, Cov, FPKM, TPM___AT1G01010, SRX1892482, 16.377369, 5.754414, 8.923347
my (%geneToFieldTypeToTpmList, $geneToFieldTypeToTpmListHref, $geneId, $experimentId, $fieldType, $tpm, $cov, $fpkm);
my ($fieldNameList, $fieldValueList, @fieldName, @fieldValue, $attrName);
# $specifiedField为关注哪一列，比如Tissue还是Treatment
# $attrName为上述指定列下的实验的类别，比如Tissue下的leaf
$geneToFieldTypeToTpmListHref = \%geneToFieldTypeToTpmList;
open FF, "<$geneExpressTsv";
while($line=<FF>){
	chomp($line);
	($fieldNameList, $fieldValueList) = split(/___/, $line);
	($geneId, $experimentId, $cov, $fpkm, $tpm) = split(/, /, $fieldValueList);
	
	if(exists($experimentInfoHref->{$experimentId}->{$specifiedField})){	
#		print "specified:" . $specifiedField . "\n";
		$attrName = $experimentInfoHref->{$experimentId}->{$specifiedField};
		if(not(exists($geneToFieldTypeToTpmListHref->{$geneId}->{$attrName}))){
			$geneToFieldTypeToTpmListHref->{$geneId}->{$attrName} = sprintf("%.4f", $tpm);
#			print join("\t", "geneId:" . $geneId, "tissue:" . $attrName, "tpmList:".$geneToFieldTypeToTpmListHref->{$geneId}->{$attrName});
		}else{
			$geneToFieldTypeToTpmListHref->{$geneId}->{$attrName} .= "," . sprintf("%.4f", $tpm);
#			print join("\t", "geneId:" . $geneId, "tissue:" . $attrName, "tpmList:".$geneToFieldTypeToTpmListHref->{$geneId}->{$attrName});
#			<STDIN>;
		}
	}else{
		next;
	}
}
close FF;


# 依次审查每个基因，挑选：
# (1) 总的组织（处理）数量达到5个以上；
# (2) 每个组织（处理）下实验样品数量达到5个以上
# (3) 每个组织（处理）下实验样本的变异系数小于0.5
# (4) 将该基因下所有组织（处理）平均表达值合在一起计算变异系数

my (@geneId, @attrName, @express, $avgExpress, $exprList);
my (%firstValidedGeneToFieldtypeToTpmList, $firstValidedGeneToFieldtypeToTpmListHref, $validedTypeNum, $cv, $outputLine, $exptNum);
$firstValidedGeneToFieldtypeToTpmListHref = \%firstValidedGeneToFieldtypeToTpmList;
@geneId = keys(%geneToFieldTypeToTpmList);

open WW, ">$outputExprCVMatrix";
print WW join("\t", "geneId", "CV", "avgExpr", "1stTissue/Treatment|exptNum|avgTpm", "2ndTissue/Treatment|exptNum|avgTpm", "3rdTissue/Treatment|exptNum|avgTpm", "...", "avgExprList", "allExprList") . "\n";
foreach $geneId(@geneId){

	# 获得该gene下收集的所有attrName，如leaf,root,flower等
	@attrName = ();
	@attrName = keys(%{$geneToFieldTypeToTpmListHref->{$geneId}});
	
####   1  第1次判断组织类型数量是否达标  ########
	#该基因下总的组织类型数量不达标(如组织类型数量没有达到5个)，那么该组织不用于识别稳定表达基因
	next if($#attrName < $minFieldTypeNum);

####   2  收集每种组织（处理）类型下的表达值并计算平均值  ######

	# 收集所有组织（处理）类型下基因表达平均值，并且计算变异系数
	# 控制变异系数不大于0.5的组织（处理）
	%firstValidedGeneToFieldtypeToTpmList=();
	$validedTypeNum = 0;
	foreach $attrName(@attrName){
		# 判断当前组织（处理）类型(比如leaf)下，实验样本数是否达标
		@express = ();
		@express = split(/,/, $geneToFieldTypeToTpmListHref->{$geneId}->{$attrName});

		# 当前组织（处理），实验样本数（表达值数）不达标，则放弃该组织
		next if($#express+1 < $minReplicateNum);

		# 检测当前组织或者处理下实验样本变异严重，放弃该组织（处理）参与后续稳定表达基因的识别
		next if(&getCV($geneToFieldTypeToTpmListHref->{$geneId}->{$attrName}) > 0.5);

		# 否则登记到%firstValidedGeneToFieldtypeToTpmist
		$validedTypeNum++;
		$firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"tpmList"} = $geneToFieldTypeToTpmListHref->{$geneId}->{$attrName};

		# 计算该组织（处理）类型下基因表达平均值和有效实验样本个数
		$avgExpress = sprintf("%.3f", sum(@express)/($#express+1));
		$firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"avgTpm"} = $avgExpress;
		$firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"tpmNum"} = $#express+1;
	}
	
	# 进一步检查该基因下有效组织（处理）类型的数量是否达标
	next if($validedTypeNum < $minFieldTypeNum);

	# 将所有有效组织（处理）类型下的所有实验样本中的表达值及平均值收集起来
	$outputLine = "";
	my $avgExprList = "";
	my $allExprList = "";

        # 列出所有有效组织（处理）名称
        @attrName = ();
        @attrName = keys(%{$firstValidedGeneToFieldtypeToTpmListHref->{$geneId}});
	foreach $attrName(@attrName){
		$allExprList .= $firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"tpmList"} . ",";
		$avgExprList .= $firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"avgTpm"} . ",";
		# 收集每个组织（处理）下的实验样本数量和平均表达值
		$exptNum = $firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"tpmNum"};
		$avgExpress = $firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"avgTpm"};
		$outputLine .= $attrName . "|" . $exptNum . "|" . $avgExpress . "\t";
	}
	$avgExprList = substr($avgExprList, 0 , length($avgExprList) - 1);
	$allExprList = substr($allExprList, 0 , length($allExprList) - 1);
	$outputLine = substr($outputLine, 0, length($outputLine) -1);


	# 用所有组织（处理）的表达平均值计算变异系数
	($avgExpress, $cv)=(-1, -1);
	&calCov($avgExprList, \$avgExpress, \$cv);
	$cv = sprintf("%.4f", $cv);
	# 输出恒定表达的基因（一直不表达的基因不输出）
	if($avgExpress>0){
        	print WW join("\t", $geneId, $cv, $avgExpress, $outputLine, $avgExprList, $allExprList) . "\n";
	}
}
close WW;


my $totalAsNum=`cat $outputExprCVMatrix |wc -l`;
my $quarterNum=int($totalAsNum/4);

my $cmd = "echo -e \"geneId\tCV\tAvgExpress\" > $outputConstant";
system($cmd);
$cmd = "grep -v \"geneId\" $outputExprCVMatrix | sort -k2,2n |awk -F \'\\t\' \'{print \$1\"\\t\"\$2\"\\t\"\$3}\' | head -n $quarterNum >> $outputConstant";
system($cmd);


# 计算一组数据的变异系数，返回这组数据平均值和变异系数值
sub calCov{
        my ($exprList,  $avgexpr, $cov) = @_;
        my (@expr, $expr, $exprU, $sum, $exprVar);
        
        @expr = split(/,/, $exprList);
        $exprU = (sum(@expr))/($#expr+1);
        $$avgexpr = $exprU;
        if($exprU == 0){
                $$cov = 0;
                return;
        }
        
        $sum = 0;
        foreach $expr(@expr){
                $sum+=($expr-$exprU)*($expr-$exprU)
        }
        $exprVar = ($sum/($#expr))**0.5;
        
        $$cov = $exprVar/$exprU;
}

sub getCV{
        my ($exprList) = @_;
        my (@expr, $expr, $exprU, $sum, $exprVar, $avgexpr, $cv);
        
        @expr = split(/,/, $exprList);
        $exprU = (sum(@expr))/($#expr+1);
        $avgexpr = $exprU;
        if($exprU == 0){
                $cv = 0;
                return $cv;
        }
        
        $sum = 0;
        foreach $expr(@expr){
                $sum+=($expr-$exprU)*($expr-$exprU)
        }
        $exprVar = ($sum/($#expr))**0.5;
        
        $cv = $exprVar/$exprU;
	return $cv;
}
