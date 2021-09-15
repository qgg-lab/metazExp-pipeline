#!/usr/bin/perl
use strict;
use Statistics::R;
use List::Util qw/max min sum maxstr minstr shuffle/;
use Getopt::Long;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--abundanceDir  \\\n" .
                "--experimentFile \\\n" .
		"--ignoreTypeList seedling,mix,unknown \\\n" .
		"--specifiedField Tissue\\\n" .
                "--minReplicateNum \\\n" .
                "--minFieldTypeNum \\\n" .
		"--foldChange 2\\\n" .
		"--outputValidedGeneExpressList \n";
	exit;
}

my ($specifiedField, $experimentFile, $geneExpressTsv, $foldChange,
$minReplicateNum, $minFieldTypeNum, $ignoreTypeList, $abundanceDir,
$outputValidedGeneExpressList);

GetOptions(
        'experimentFile=s'=>\$experimentFile,
	'ignoreTypeList=s'=>\$ignoreTypeList,
	'specifiedField=s'=>\$specifiedField,
        'abundanceDir=s'=>\$abundanceDir,
	'foldChange=s'=>\$foldChange,
        'minReplicateNum=s'=>\$minReplicateNum,
        'minFieldTypeNum=s'=>\$minFieldTypeNum,
	'outputValidedGeneExpressList=s'=>\$outputValidedGeneExpressList,
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

	# 如果该组织/处理是关注的类型，那么放弃登记
	next if(index($ignoreTypeList, $tmpHash{$specifiedField})>=0);

	$experimentInfoHref->{$tmpHash{"Experiment"}}->{$specifiedField} = $tmpHash{$specifiedField};
}
close FF;

## == 2 构建gene->fieldType->tpmList ===
# GeneId, ExptId, Cov, FPKM, TPM___AT1G01010, SRX1892482, 16.377369, 5.754414, 8.923347
my (%geneToFieldTypeToTpmList, $geneToFieldTypeToTpmListHref, $geneId, @experimentId, $experimentId, $fieldType, $tpm, $cov, $fpkm,  $abundanceFile);
my ($fieldNameList, $fieldValueList, @fieldName, @fieldValue, $attrName);
# $specifiedField为关注哪一列，比如Tissue还是Treatment
# $attrName为上述指定列下的实验的类别，比如Tissue下的leaf
$geneToFieldTypeToTpmListHref = \%geneToFieldTypeToTpmList;

@experimentId = keys(%exprimentInfo);
foreach $experimentId(@experimentId){
	$abundanceFile = $abundanceDir . "/" . $experimentId . "/geneAbundanceByStringtie.tab";
	# Gene ID Gene Name       Reference       Strand  Start   End     Coverage        FPKM    TPM
	# AT1G01010       NAC001  1       +       3631    5899    6.687401        3.133270        4.926901
	open FF, "<$abundanceFile";
	<FF>;
	while($line=<FF>){
		chomp($line);
		@fieldValue = ();
		@fieldValue = split(/\t/, $line);
		($geneId, $cov, $fpkm, $tpm) = ($fieldValue[0], $fieldValue[6], $fieldValue[7], $fieldValue[8]);
		$attrName = $experimentInfoHref->{$experimentId}->{$specifiedField};
		# 将TPM值添加到基因下对应的组织类型中
		# 如：$geneToFieldTypeToTpmListHref->{"AT1G16010"}->{"leaf"}="123,90,45,23,0.3";
		if(not(exists($geneToFieldTypeToTpmListHref->{$geneId}->{$attrName}))){
			$geneToFieldTypeToTpmListHref->{$geneId}->{$attrName} = sprintf("%.4f", $tpm);
		}else{
			$geneToFieldTypeToTpmListHref->{$geneId}->{$attrName} .= "," . sprintf("%.4f", $tpm);
		}
	}
	close FF;
}

# 依次审查每个基因，挑选：
# (1) 总的实验类型（处理）数量达到5个以上；
# (2) 每个实验类型（处理）下重复样品数量达到5个以上
# (3) 找到表达水平平均值最高的实验类型
# (4) 找到表达水平平均值最低的实验类型
# (5) 将符合上述两个条件的基因及其在各种类型实验中的表达值输出到outputValidedGeneExpressList
# specTag0.01  specTag0.05   geneId direction   specialType   specialAvgExpress	tissueExprPvalue    pvalueTissue2 pvalueTissue3
#	Y      Y             ATG102  High	leaf	      357.69		root|112|0.002	 stem|80|0.005    seed|34|0.567
# ===

# $specifiedField=Tissue/Treatment
# $attrName=leaf

# 建立一个R连接
my ($R, $cmd);
$R = Statistics::R->new( shared => 1);
my (@geneId, @attrName, @express, $avgExpress);
my (%firstValidedGeneToFieldtypeToTpmList, $firstValidedGeneToFieldtypeToTpmListHref, $validedTypeNum, $maxExpressAttrName, $minExpressAttrName, $maxExpress, $minExpress);
$firstValidedGeneToFieldtypeToTpmListHref = \%firstValidedGeneToFieldtypeToTpmList;
@geneId = keys(%geneToFieldTypeToTpmList);
open WW, ">$outputValidedGeneExpressList";
print WW join("\t", "pvalue0.01Rlt", "pvalue0.03Rlt", "pvalue0.05Rlt", "geneId", "regulation", "tissue/treatment", "avgExpTpm", "experimentNum", "1stTissue/treatment|avgTpm|exptNum|pvalue", "2ndTissue/treatment|avgTpm|exptNum|pvalue", "3rdTissue/Treatment|avgTpm|exptNum|pvalue", "4thTissue/Treatment|avgTpm|exptNum|pvalue", "...") . "\n";
foreach $geneId(@geneId){

	# 获得该gene下收集的所有attrName，如leaf,root,flower等
	# $geneToFieldTypeToTpmListHref->{"AT1G16010"} => @attrName = ("leaf", "root", "meristem", "flower", "seed")
	@attrName = ();
	@attrName = keys(%{$geneToFieldTypeToTpmListHref->{$geneId}});
	
####   1  第1次判断实验类型数量是否达标  ########
	#该基因下总的实验类型数量不达标(如组织类型数量没有达到5个)
	next if($#attrName+1 < $minFieldTypeNum);

####   2  收集每种实验类型下的表达值并计算平均值  ######
	# 收集所有实验类型下基因表达值列表
	%firstValidedGeneToFieldtypeToTpmList=();
	$validedTypeNum = 0;
	foreach $attrName(@attrName){
		# 判断当前实验类型(比如leaf)下，实验重复数是否达标
		@express = ();
		@express = split(/,/, $geneToFieldTypeToTpmListHref->{$geneId}->{$attrName});
		# $geneToFieldTypeToTpmListHref->{"AT1G16010"}->{"leaf"} = "123,90,45,23,0.3"

		# 当前实验类型(如leaf)下，实验重复数（表达值数）不达标，则放弃该实验类型
		next if($#express+1<$minReplicateNum);

		# 否则登记到%firstValidedGeneToFieldtypeToTpmist
		# 对不符合重复数要求的组织类型/处理类型放弃后，将剩余的组织类型/处理类型的表达值转移到firstValiddGeneTo……中
		# $firstValidedGeneToFieldtypeToTpmListHref->{"AT1G16010"}->{"leaf"}->{"tpmList"} = "123,90,45,23,0.3";
		$validedTypeNum++;
		$firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"tpmList"} = $geneToFieldTypeToTpmListHref->{$geneId}->{$attrName};

		# 计算该实验类型下基因表达平均值
		# $firstValidedGeneToFieldtypeToTpmListHref->{"AT1G16010"}->{"leaf"}->{"avgTpm"} = 87.23
		# $firstValidedGeneToFieldtypeToTpmListHref->{"AT1G16010"}->{"leaf"}->{"tpmNum"} = 5;
		$avgExpress = sprintf("%.3f", sum(@express)/($#express+1));
		$firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"avgTpm"} = $avgExpress;
		$firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"tpmNum"} = $#express+1;
	}
	
	# 进一步检查该基因下符合要求的指定类型数量是否达标
	next if($validedTypeNum < $minFieldTypeNum);




######### 准备找表达平均值最高和最低的组织/处理  #################################
	$maxExpressAttrName = "";
	$maxExpress = -1;
	$minExpressAttrName = "";
	$minExpress = 100000;

	# 遍历所有attrName找到最高和最低表达类型及平均表达值
	@attrName = ();
	@attrName = keys(%{$firstValidedGeneToFieldtypeToTpmListHref->{$geneId}});
	foreach $attrName(@attrName){
		if($maxExpress < $firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"avgTpm"}){
			$maxExpress = $firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"avgTpm"};
			$maxExpressAttrName = $attrName;
		}
		if($minExpress > $firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"avgTpm"}){
			$minExpress = $firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"avgTpm"};
			$minExpressAttrName = $attrName;
		}
	}

###################################################################################################
#												  #
#                               输出最高和最低表达的组织/处理
#                               并计算最高和最低表达的组织/处理是否显著                           #
#                               								  #
###################################################################################################
	my ($outputHighLine, $outputLowLine) =  ("", "");
	my ($expressList1, $expressList2, $pvalue, $pvalue001Flag, $pvalue003Flag, $pvalue005Flag, $avgExpress, $tpmNum);

	#################################################################################
	#										#
	#				最高表达组织/处理输出				#
	#										#
	#################################################################################

	# 获得最高表达组织/处理上Tpm表达值列表
	$expressList1 = $firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$maxExpressAttrName}->{"tpmList"};
	$tpmNum = $firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$maxExpressAttrName}->{"tpmNum"};
	my $avgTpm1 =$firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$maxExpressAttrName}->{"avgTpm"};
	# 输出最高表达的组织/处理的名称和平均表达水平
	($pvalue001Flag, $pvalue003Flag, $pvalue005Flag) = ("Y", "Y", "Y");
	$outputHighLine = join("\t", $geneId, "High", $maxExpressAttrName, $maxExpress, $tpmNum);
	foreach $attrName(@attrName){
		next if($attrName eq $maxExpressAttrName);
		# 获得1个非最高表达组织上TPM表达值列表，和平均表达值
		$expressList2 = $firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"tpmList"};
		$avgExpress = $firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"avgTpm"};
		$tpmNum = $firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"tpmNum"};

		# 计算最高表达组织/处理和当前组织/处理的TPM之间的差异显著水平，放在pvalue中
$cmd=<<EOF;
x<-c($expressList1)
y<-c($expressList2)
pvalue<-t.test(x, y, alternative = "greater", paired=FALSE)\$p.value
EOF
#print "R cmd:" . $cmd;
#<SDTIN>;


		if(&innerSame($expressList1)==1 and &innerSame($expressList2)==1 and $avgTpm1 == $avgExpress){
		# 两组数据内部都相同，组间平均值也相同
			$pvalue = 1;
		}elsif(&innerSame($expressList1)==1 and &innerSame($expressList2)==1 and $avgTpm1 != $avgExpress){
		# 两组数据内部都相同，组间平均值不相同
			$pvalue = 0;
		}elsif(&innerSame($expressList1)!=1 or &innerSame($expressList2)!=1){
		# 两组数据至少1组内部不相同
			$R->run($cmd);
			$pvalue=$R->get('pvalue');
		}
		# 登记是否真的特异高，如果有1个pvalue大于001/0.05，那么就不显著
		$pvalue001Flag = "N" if($pvalue > 0.01);
		$pvalue003Flag = "N" if($pvalue > 0.03);
		$pvalue005Flag = "N" if($pvalue > 0.05);

		# 检测倍数
		if($avgExpress!=0 and $avgTpm1/$avgExpress < $foldChange){
			($pvalue001Flag, $pvalue003Flag, $pvalue005Flag) = ("N", "N", "N");
		}elsif($avgExpress==0 and $avgTpm1==0){
			($pvalue001Flag, $pvalue003Flag, $pvalue005Flag) = ("N", "N", "N");
		}

		# 将当前非最高表达组织上的avgTPM和pvalue输出
		$outputHighLine .= "\t" . $attrName . "|" . $avgExpress . "|" . $tpmNum .  "|" . sprintf("%.3f", $pvalue);
	}


	# 输出最高表达组织和其它组织的显著差异水平
	print WW join("\t", $pvalue001Flag, $pvalue003Flag, $pvalue005Flag, $outputHighLine) . "\n";




	#################################################################################
	#										#
	#				最低表达组织/处理输出				#
	#										#
	#################################################################################

	# 获得最低表达组织/处理上Tpm表达值列表
	$expressList1 = $firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$minExpressAttrName}->{"tpmList"};
	$tpmNum = $firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$minExpressAttrName}->{"tpmNum"};
	my $avgTpm1 =$firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$minExpressAttrName}->{"avgTpm"};

	# 输出最低表达组织/处理上的名称和平均表达水平
	($pvalue001Flag, $pvalue003Flag, $pvalue005Flag) = ("Y", "Y", "Y");
	$outputLowLine = join("\t", $geneId, "Low", $minExpressAttrName, $minExpress, $tpmNum);
	foreach $attrName(@attrName){
		next if($attrName eq $minExpressAttrName);

		# 获得1个最低表达组织/处理上TPM表达值列表，和平均表达值
		$expressList2 = $firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"tpmList"};
		$avgExpress = $firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"avgTpm"};
		$tpmNum = $firstValidedGeneToFieldtypeToTpmListHref->{$geneId}->{$attrName}->{"tpmNum"};	

		# 计算最低表达组织/处理和当前的TPM之间的差异显著水平，放在pvalue中
$cmd=<<EOF;
x<-c($expressList1)
y<-c($expressList2)
pvalue<-t.test(x, y, alternative = "less", paired=FALSE)\$p.value
EOF

		if(&innerSame($expressList1)==1 and &innerSame($expressList2)==1 and $avgTpm1 == $avgExpress){
		# 两组数据内部都相同，组间平均值也相同
			$pvalue = 1;
		}elsif(&innerSame($expressList1)==1 and &innerSame($expressList2)==1 and $avgTpm1 != $avgExpress){
		# 两组数据内部都相同，组间平均值不相同
			$pvalue = 0;
		}elsif(&innerSame($expressList1)!=1 or &innerSame($expressList2)!=1){
		# 两组数据至少1组内部不相同
			$R->run($cmd);
			$pvalue=$R->get('pvalue');
		}

		# 登记是否真的特异低，如果有1个pvalue大于001/003/005，那么就不显著
		$pvalue001Flag = "N" if($pvalue > 0.01);
		$pvalue003Flag = "N" if($pvalue > 0.03);
		$pvalue005Flag = "N" if($pvalue > 0.05);

		# 检测倍数
		if($avgTpm1!=0 and $avgExpress/$avgTpm1 < $foldChange){
			($pvalue001Flag, $pvalue003Flag, $pvalue005Flag) = ("N", "N", "N");
		}elsif($avgExpress==0 and $avgTpm1==0){
			($pvalue001Flag, $pvalue003Flag, $pvalue005Flag) = ("N", "N", "N");
		}

		# 将当前非特异低表达组织上的TPM添加的输出行: root|345|0.04
		$outputLowLine .= "\t" . $attrName . "|" . $avgExpress . "|" . $tpmNum . "|" . sprintf("%.5f", $pvalue);
	}

	# 输出最地表达组织和其它组织的显著差异水平
	print WW join("\t", $pvalue001Flag, $pvalue003Flag, $pvalue005Flag, $outputLowLine) . "\n";	
	
}
close WW;


$R->stop;


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

