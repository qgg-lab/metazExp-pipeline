#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--python \\\n" .
                "--rmatsStat \\\n" .
                "--asReadCountDir \\\n" .
                "--controlExptIdList \\\n" .
                "--treatmentExptIdList \\\n" .
                "--diff_cutoff \\\n" .
                "--threadNum 4 \\\n" .
		"--outputDir \n";
	exit;
}

my ($python, $rmatsStat, $asReadCountDir, $controlExptIdList, $treatmentExptIdList, $diff_cutoff, $threadNum, $outputDir);

GetOptions(
        'python=s'=>\$python,
        'rmatsStat=s'=>\$rmatsStat,
        'asReadCountDir=s'=>\$asReadCountDir,
        'controlExptIdList=s'=>\$controlExptIdList,
	'treatmentExptIdList=s'=>\$treatmentExptIdList,
	'diff_cutoff=s'=>\$diff_cutoff,
	'threadNum=s'=>\$threadNum,
	'outputDir=s'=>\$outputDir,
);

my (@controlExptId, @treatmentExptId, $exptId);
my (%asReadCountInfo, $asReadCountInfoHref);
$asReadCountInfoHref = \%asReadCountInfo;


# 输入：
# ASID                  IJC_SAMPLE_1    SJC_SAMPLE_1    IncFormLen      SkipFormLen     Experiment
# OSJGA3SS0000000001    0               1               78              74              ERX2120647

# 输出：
# ASID			IJC1		SJC1		IJC2		SJC2		IncFormLen	SkipFormLen
# OSJGA3SS0000000001	0,1,3		2,3,4		1,5,9		9,3,4		78		74

# 扫描每一个expt的read情况，登记到as->expt->{IJC/SJC}/
my (@field, %asId, @asId, %IncFormLen, %SkipFormLen, @IncFormLen, @SkipFormLen);
@controlExptId = ();
@controlExptId = split(/,/, $controlExptIdList);
foreach $exptId(@controlExptId){
	open FF, "<$asReadCountDir/$exptId";
	<FF>;
	while(my $line=<FF>){
		chomp($line);
		@field = split(/\t/, $line);
		# 放弃没有检测到read覆盖的AS事件
		if($field[1]!=0 or $field[2]!=0){
			$asId{$field[0]} = 1;
			$asReadCountInfoHref->{"control"}->{$exptId}->{$field[0]}->{"IJC"} = $field[1];
			$asReadCountInfoHref->{"control"}->{$exptId}->{$field[0]}->{"SJC"} = $field[2];
			$asReadCountInfoHref->{"control"}->{$exptId}->{$field[0]}->{"IncFormLen"} = $field[3];
			$asReadCountInfoHref->{"control"}->{$exptId}->{$field[0]}->{"SkipFormLen"} = $field[4];
		}
	}
	close FF;
}

@treatmentExptId = ();
@treatmentExptId = split(/,/, $treatmentExptIdList);
foreach $exptId(@treatmentExptId){
	open FF, "<$asReadCountDir/$exptId";
	<FF>;
	while(my $line=<FF>){
		chomp($line);
		@field = split(/\t/, $line);
		# 放弃没有检测到read覆盖的AS事件
		if($field[1]!=0 or $field[2]!=0){
			$asId{$field[0]} = 1;
			$asReadCountInfoHref->{"treatment"}->{$exptId}->{$field[0]}->{"IJC"} = $field[1];
			$asReadCountInfoHref->{"treatment"}->{$exptId}->{$field[0]}->{"SJC"} = $field[2];
			$asReadCountInfoHref->{"treatment"}->{$exptId}->{$field[0]}->{"IncFormLen"} = $field[3];
			$asReadCountInfoHref->{"treatment"}->{$exptId}->{$field[0]}->{"SkipFormLen"} = $field[4];
		}
	}
	close FF;
}


# 分别提取control和treatment的expt中的read count，然后联立起来输出
my ($outputLine, $asId);
my ($allControlExits, $allTreatmentExits, $controlIjcList, $controlSjcList, $treatmentIjcList, $treatmentSjcList);
@asId = keys(%asId);
system("mkdir -p $outputDir");
my $incTxt = $outputDir . "/inc.txt";
open WW, ">$incTxt";
print WW join("\t", "Exon", "IJC1", "SJC1", "IJC2", "SJC2", "IncFormLen", "SkipFormLen") . "\n";
foreach $asId(@asId){

	$outputLine = $asId;
	
	%IncFormLen = ();
	%SkipFormLen = ();
	@IncFormLen = ();
	@SkipFormLen = ();

	# 提取control中IJC和ISC
	$allControlExits = 1;
	$controlIjcList = "";
	$controlSjcList = "";
	foreach $exptId(@controlExptId){
		if(exists($asReadCountInfoHref->{"control"}->{$exptId}->{$asId}->{"IJC"})){
			$controlIjcList .= $asReadCountInfoHref->{"control"}->{$exptId}->{$asId}->{"IJC"} . ",";
			$IncFormLen{$asReadCountInfoHref->{"control"}->{$exptId}->{$asId}->{"IncFormLen"}} = 1;
			$SkipFormLen{$asReadCountInfoHref->{"control"}->{$exptId}->{$asId}->{"SkipFormLen"}} = 1;
		}else{
			$allControlExits = 0;
		}

		if(exists($asReadCountInfoHref->{"control"}->{$exptId}->{$asId}->{"SJC"})){
			$controlSjcList .= $asReadCountInfoHref->{"control"}->{$exptId}->{$asId}->{"SJC"} . ",";
			$IncFormLen{$asReadCountInfoHref->{"control"}->{$exptId}->{$asId}->{"IncFormLen"}} = 1;
			$SkipFormLen{$asReadCountInfoHref->{"control"}->{$exptId}->{$asId}->{"SkipFormLen"}} = 1;
		}else{
			$allControlExits = 0;
		}
	}
	
	# 提取treatment中的IJC和ISC
	$allTreatmentExits = 1;
 	$treatmentIjcList = "";
	$treatmentSjcList = "";	
	foreach $exptId(@treatmentExptId){
		if(exists($asReadCountInfoHref->{"treatment"}->{$exptId}->{$asId}->{"IJC"})){
			$treatmentIjcList .= $asReadCountInfoHref->{"treatment"}->{$exptId}->{$asId}->{"IJC"} . ",";
			$IncFormLen{$asReadCountInfoHref->{"treatment"}->{$exptId}->{$asId}->{"IncFormLen"}} = 1;
			$SkipFormLen{$asReadCountInfoHref->{"treatment"}->{$exptId}->{$asId}->{"SkipFormLen"}} = 1;
		}else{
			$allControlExits = 0;
		}
		if(exists($asReadCountInfoHref->{"treatment"}->{$exptId}->{$asId}->{"SJC"})){
			$treatmentSjcList .= $asReadCountInfoHref->{"treatment"}->{$exptId}->{$asId}->{"SJC"} . ",";
			$IncFormLen{$asReadCountInfoHref->{"treatment"}->{$exptId}->{$asId}->{"IncFormLen"}} = 1;
			$SkipFormLen{$asReadCountInfoHref->{"treatment"}->{$exptId}->{$asId}->{"SkipFormLen"}} = 1;
		}else{
			$allTreatmentExits = 0;
		}	
	}

	# 检查该AS下所有experimnt中采集的IncFormLen和SkipFormLen是否唯一
	@IncFormLen = keys(%IncFormLen);
	@SkipFormLen = keys(%SkipFormLen);
	next if($#IncFormLen>0 or $#SkipFormLen>0);
	
	# 检查该AS是否在所有的control和treatment中都被read覆盖
	# 该要求是否强制要求？需要进一步检测
	next if($allControlExits==0 or $allTreatmentExits==0);

	# 输出
	print WW join("\t", $asId, substr($controlIjcList, 0, length($controlIjcList) - 1), substr($controlSjcList, 0, length($controlSjcList) - 1), substr($treatmentIjcList, 0, length($treatmentIjcList) - 1), substr($treatmentSjcList, 0, length($treatmentSjcList) - 1), $IncFormLen[0], $SkipFormLen[0]) . "\n";

}
close WW;

# 执行差异分析
my $cmd = "$python $rmatsStat $incTxt $outputDir $threadNum $diff_cutoff &>/dev/null";
system($cmd);
