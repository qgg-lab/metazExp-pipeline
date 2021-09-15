#!/usr/bin/perl
use strict;
use List::Util qw/max min sum maxstr minstr shuffle/;
use Getopt::Long;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--experimentFile \\\n" .
		"--unitVolum \\\n" .
		"--minReadCovPerUnitVolum \\\n" .
                "--psiFileList \\\n" .
                "--minReplicateNum \\\n" .
                "--minTissueNum \\\n" .
                "--minDltPsi \\\n" .
                "--outputTissueSpecificHighPsiAs \\\n" .
		"--outputTissueSpecificLowPsiAs \n";
	exit;
}

my ($experimentFile, $unitVolum, $minReadCovPerUnitVolum, $psiFileList, 
$minReplicateNum, $minTissueNum, $minDltPsi, $outputTissueSpecificHighPsiAs, $outputTissueSpecificLowPsiAs);

GetOptions(
        'experimentFile=s'=>\$experimentFile,
        'unitVolum=s'=>\$unitVolum,
        'minReadCovPerUnitVolum=s'=>\$minReadCovPerUnitVolum,
        'psiFileList=s'=>\$psiFileList,
        'minReplicateNum=s'=>\$minReplicateNum,
        'minTissueNum=s'=>\$minTissueNum,
	'minDltPsi=s'=>\$minDltPsi,
	'outputTissueSpecificHighPsiAs=s'=>\$outputTissueSpecificHighPsiAs,
	'outputTissueSpecificLowPsiAs=s'=>\$outputTissueSpecificLowPsiAs,
);
my ($line, $i, @nameField, @valueField);
##  == 1 获得experiment的组织信息 ===
# 将experiment读入hash，可以通过experimentId直接获得tissue类型
my (%exprimentInfo, $experimentInfoHref, %tmpHash, $exptId, $tissue);
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
	$experimentInfoHref->{$tmpHash{"Experiment"}}->{"tissue"} = $tmpHash{"Tissue"};
	$experimentInfoHref->{$tmpHash{"Experiment"}}->{"minReadCoverage"} = $tmpHash{"mappedBases"}/$unitVolum * $minReadCovPerUnitVolum;
}
close FF;


## == 2 构建as->tissue->psiList hash ===
# ASID    IJC_SAMPLE_1    SJC_SAMPLE_1    IncFormLen      SkipFormLen     Experiment
# ATHAA3SS0000000001      549     189     398     100     SRX485073
my (%asToTissueToPsi, $asToTissueToPsiHref, $asId, $experimentId, $tissue);
my ($inclusionNorm, $exclusionNorm, $psi);
my (@psiFile, $psiFile);
$asToTissueToPsiHref = \%asToTissueToPsi;
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
		$tissue = $experimentInfoHref->{$experimentId}->{"tissue"};
		$asId = $tmpHash{"ASID"};

		$inclusionNorm = $tmpHash{"IJC_SAMPLE_1"}/$tmpHash{"IncFormLen"};
		$exclusionNorm = $tmpHash{"SJC_SAMPLE_1"}/$tmpHash{"SkipFormLen"};
		$psi=$inclusionNorm/($inclusionNorm+$exclusionNorm);

		# 将Psi值登记到AS的相关tissue中
		if(not(exists($asToTissueToPsiHref->{$asId}->{$tissue}))){
			$asToTissueToPsiHref->{$asId}->{$tissue} = $psi;
		}else{
			$asToTissueToPsiHref->{$asId}->{$tissue} .= "," . $psi;
		}
	}
	close FF;
}

# 输出包含组织特异的as
my (@asId, @tissue, $checkReplicate, @psi, $validTissueNum, 
%validTissueToAvgPsi, $validTissueToAvgPsiHref, 
$maxTissue, $minTissue, $maxPsi, $minPsi, 
$checkSpecificMax, $checkSpecificMin);
open MAXWW, ">$outputTissueSpecificHighPsiAs";
open MINWW, ">$outputTissueSpecificLowPsiAs";
print MAXWW join("\t", "ASID", "Tissue=MaxPsi", "AllTissuePsiList") . "\n";
print MINWW join("\t", "ASID", "Tissue=MinPsi", "AllTissuePsiList") . "\n";
@asId = keys(%asToTissueToPsi);
foreach $asId(@asId){
	@tissue = ();
	@tissue = keys(%{$asToTissueToPsiHref->{$asId}});
	# 检查每个tissue中收集的PSI数量是否达到指定值: minReplicateNum
	$validTissueNum = 0;
	%validTissueToAvgPsi = ();
	$validTissueToAvgPsiHref = \%validTissueToAvgPsi;
	foreach $tissue(@tissue){
		@psi = ();
		@psi = split(",", $asToTissueToPsiHref->{$asId}->{$tissue});
		if($#psi+1 >= $minReplicateNum){
			$validTissueNum++;
			# 计算该组织内的psi平均值
			$validTissueToAvgPsiHref->{$tissue} = sum(@psi)/($#psi+1);
		}
	}

	# 检查有效tissue数量是否达到指定值
	next if($validTissueNum < $minTissueNum);

	# 找出最大psi平均值和最小psi平均值对应的组织
	$maxPsi = -1;
	$minPsi = 2;
	@tissue = ();
	@tissue = keys(%validTissueToAvgPsi);
	foreach $tissue(@tissue){
		if($maxPsi < $validTissueToAvgPsiHref->{$tissue}){
			$maxPsi = $validTissueToAvgPsiHref->{$tissue};
			$maxTissue = $tissue;
		}
		if($minPsi > $validTissueToAvgPsiHref->{$tissue}){
			$minPsi = $validTissueToAvgPsiHref->{$tissue};
			$minTissue = $tissue;
		}
	}

	# 检查最大psi平均值是否比其他所有组织的psi值都大0.25
	# 检查最小psi平均值是否比其他所有组织的psi值都小0.25
	$checkSpecificMax = 1;
	$checkSpecificMin = 1;
	foreach $tissue(@tissue){
		if($maxTissue ne $tissue and $validTissueToAvgPsiHref->{$maxTissue} - $validTissueToAvgPsiHref->{$tissue} < $minDltPsi){
			$checkSpecificMax = 0;
		}
		if($minTissue ne $tissue and $validTissueToAvgPsiHref->{$tissue} - $validTissueToAvgPsiHref->{$minTissue} < $minDltPsi){
			$checkSpecificMin = 0;
		}
	}

	# 如果maxTissue和minTissue检查都没通过，那么放弃
	if($checkSpecificMax == 0 and $checkSpecificMin == 0){
		next;
	}elsif($checkSpecificMax == 1 and $checkSpecificMin == 1){
		print MAXWW join("\t", $asId, $maxTissue . "=" . sprintf("%.3f", $validTissueToAvgPsiHref->{$maxTissue}));
		foreach $tissue(@tissue){
			if($tissue ne $maxTissue){
				print MAXWW "\t" . $tissue . "=" . sprintf("%.3f", $validTissueToAvgPsiHref->{$tissue});
			}
		}
		print MAXWW "\n";

		print MINWW join("\t", $asId, $minTissue . "=" . sprintf("%.3f", $validTissueToAvgPsiHref->{$minTissue}));
		foreach $tissue(@tissue){
			if($tissue ne $minTissue){
				print MINWW "\t" . $tissue . "=" . sprintf("%.3f", $validTissueToAvgPsiHref->{$tissue});
			}
		}
		print MINWW "\n";
	}elsif($checkSpecificMax == 1 and $checkSpecificMin == 0){
		print MAXWW join("\t", $asId, $maxTissue . "=" . sprintf("%.3f", $validTissueToAvgPsiHref->{$maxTissue}));
		foreach $tissue(@tissue){
			if($tissue ne $maxTissue){
				print MAXWW "\t" . $tissue . "=" . sprintf("%.3f", $validTissueToAvgPsiHref->{$tissue});
			}
		}
		print MAXWW "\n";
	}elsif($checkSpecificMax == 0 and $checkSpecificMin == 1){
		print MINWW join("\t", $asId, $minTissue . "=" . sprintf("%.3f", $validTissueToAvgPsiHref->{$minTissue}));
		foreach $tissue(@tissue){
			if($tissue ne $minTissue){
				print MINWW "\t" . $tissue . "=" . sprintf("%.3f", $validTissueToAvgPsiHref->{$tissue});
			}
		}
		print MINWW "\n";

	}
}
close MAXWW;
close MINWW;
