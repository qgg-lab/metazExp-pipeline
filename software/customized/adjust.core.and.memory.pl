#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--sampleInfoTsv \\\n" .
                "--sampleIdList \\\n" .
                "--sbatchDir \n";
	exit;
}

my ($sampleInfoTsv, $sampleIdList, $sbatchDir);

GetOptions(
        'sampleInfoTsv=s'=>\$sampleInfoTsv,
	'sampleIdList=s'=>\$sampleIdList,
        'sbatchDir=s'=>\$sbatchDir,
);

my (%selectedSample, @exptId, $exptId);
open FF, "<$sampleIdList";
@exptId = <FF>;
close FF;
foreach $exptId(@exptId){
	chomp($exptId);
	$selectedSample{$exptId}=1;
}


my (%sampleVolume);
my (%exptId);
open FF, "<$sampleInfoTsv";
my (@fieldName, @fieldValue, $line, %tmpHref);
my $line = <FF>;
chomp($line);
@fieldName = split(/\t/, $line);
while($line=<FF>){
	@fieldValue = ();
	@fieldValue = split(/\t/, $line);
	%tmpHref = ();
	for(my $i=0; $i<=$#fieldValue; $i++){
		$tmpHref{$fieldName[$i]} = $fieldValue[$i];
	}
	$sampleVolume{$tmpHref{"Experiment"}}=$tmpHref{"Base"};
	$exptId{$tmpHref{"Experiment"}} = 1;
}
close FF;

my ($coreNum, $memSize, $time, $sizeBycore, $sizeByVolume, $genomeIndexSize, $enlargeVolumeCoef);
$genomeIndexSize = 2;
$enlargeVolumeCoef = 2;

@exptId = keys(%exptId);

# 每增加1G数据量，那么核数增加1
# 内存：设置核心数/28*115  ? 测序数据量G * 2 + 2G( for genome index size)
# $enlargeVolumeCoef = 2; $genomeIndexSize = 2

foreach $exptId(@exptId){

	# 需要的核数为测序容量的2倍，但是至少2核
	$coreNum = int($sampleVolume{$exptId} * 2) + 2;
	# 核数调整成偶数
	if(int($coreNum/2)*2!=$coreNum){
		$coreNum +=1;
	}

	# 需要的内存容量
	# ----按照28核计算容量推算
	$sizeBycore = int($coreNum / 28 * 115 + 8);

	# ----按照实际测序容量推算:
	$sizeByVolume = int($sampleVolume{$exptId} * $enlargeVolumeCoef + $genomeIndexSize + 8);

	# 上述取大值
	$memSize = $sizeByVolume;
	# 以28机器为参照，按照核数推算的内存 增长速度
	if($sizeBycore > $sizeByVolume){
		$memSize = $sizeBycore;
	}	
	# 把内存调整成偶数
	if(int($memSize/2)*2!=$memSize){
		$memSize +=1;
	}

	# 控制28核内的任务需要的内存不得大于115G
	if($coreNum <= 28 and $memSize>115){
		$memSize = 110;
	}

	if($memSize<32){
		$memSize = 32;
	}
	$time = "03:58:00";

	# 调整超大任务
	if($coreNum > 28 and $coreNum <40){
		$coreNum = "28";
		$memSize = "110";
		$time = "08:00:00";
	}elsif($coreNum >= 40 and $coreNum<=60){
		$coreNum = "40";
		$memSize = "178";
		$time ="12:00:00";
	}elsif($coreNum >60 and $coreNum<=80){
		$coreNum = "40";
		$memSize = "367";
		$time = "24:00:00";
	}elsif($coreNum >80){
		$coreNum = "40";
		$memSize = "367";
		$time = "36:00:00";
	}

	$memSize = $memSize . "G";

#$time = "08:00:00";	
#$memSize = "110G";
#$coreNum = "40";

	my $sbatchFile = $sbatchDir . "/" . $exptId;

	# sbatch文件不存在或者不是关注的expt
	next if(not -e $sbatchFile or not(exists($selectedSample{$exptId})));

#	print join("\t", "exptId", $exptId,  "Volume", $sampleVolume{$exptId}, "Core", $coreNum, "Memory", $memSize, "time: ", $time) . "\n";
#	print $coreNum . "\n";
#	<STDIN>;
#	next;


	my $cmd = "sed -i \'s/^\\(.*task=\\).*/\\1$coreNum/\' $sbatchFile";
	system($cmd);
	$cmd = "sed -i \'s/^\\(.*mem=\\).*/\\1$memSize/\' $sbatchFile";
	system($cmd);
	$cmd = "sed -i \'s/^\\(.*threadNum=\\).*/\\1$coreNum/\' $sbatchFile";
	system($cmd);
	$cmd = "sed -i \'s/^\\(.*threadNum \\).*\\( 1.*\\)\$/\\1$coreNum\\2/\' $sbatchFile";
	system($cmd);
	$cmd = "sed -i \'s/^\\(.*time=\\).*/\\1$time/\' $sbatchFile";
	system($cmd);

}
