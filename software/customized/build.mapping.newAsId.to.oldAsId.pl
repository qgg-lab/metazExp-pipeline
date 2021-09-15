#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--oldAsFileList A3ss,A5ss,Mxe,Ri,Se\\\n" .
                "--newAsFileList A3ss,A5ss,Mxe,Ri,Se\\\n" .
		"--jcecFileList A3ss,A5ss,Mxe,Ri,Se\\\n" .
		"--jcFileList A3ss,A5ss,Mxe,Ri,Se\\\n" .
		"--newAddedA3ssFile \\\n" .
		"--newAddedA5ssFile \\\n" .
		"--newAddedSeFile \\\n" .
		"--newAddedRiFile \\\n" .
		"--newAddedMxeFile \\\n" .
		"--abbr  \\\n" .
                "--newAsIdToUniformAsIdTsv \n";
	exit;
}

my ($oldAsFileList, $newAsFileList, $newAddedA3ssFile, $newAddedA5ssFile, $newAddedSeFile, $newAddedRiFile, $newAddedMxeFile, $newAsIdToUniformAsIdTsv, $abbr, $jcecFileList, $jcFileList);

GetOptions(
        'oldAsFileList=s'=>\$oldAsFileList,
        'newAsFileList=s'=>\$newAsFileList,
	'jcecFileList=s'=>\$jcecFileList,
	'jcFileList=s'=>\$jcFileList,
	'abbr=s'=>\$abbr,
        'newAddedA3ssFile=s'=>\$newAddedA3ssFile,
        'newAddedA5ssFile=s'=>\$newAddedA5ssFile,
        'newAddedSeFile=s'=>\$newAddedSeFile,
        'newAddedMxeFile=s'=>\$newAddedMxeFile,
        'newAddedRiFile=s'=>\$newAddedRiFile,
        'newAsIdToUniformAsIdTsv=s'=>\$newAsIdToUniformAsIdTsv,
);


# 从坐标到老AS事件ID的hash
my (%toUniformAsId);
my (@oldAsFile, $oldAsFile, @newAsFile, $newAsFile, $line, @tmp, $asId, $coordinateString);
my (%coordinateToAsId, $coordinateToAsIdHref);
$coordinateToAsIdHref = \%coordinateToAsId;
@oldAsFile = split(/,/, $oldAsFileList);
foreach $oldAsFile(@oldAsFile){
	# print $oldAsFile . "\n";
	
	open FF, "<$oldAsFile";
	<FF>;
	while($line=<FF>){
		chomp($line);
		@tmp = ();
		@tmp = split(/\t/, $line);
		$asId = shift(@tmp);
		$coordinateString = join("\t", @tmp);
		$coordinateToAsIdHref->{$coordinateString}->{"oldAsId"} = $asId;
	}
	close FF;
}


# 从坐标到新AS事件ID的hash
@newAsFile = split(/,/, $newAsFileList);
foreach $newAsFile(@newAsFile){
	#print $newAsFile . "\n";
	open FF, "<$newAsFile";
	<FF>;
	while($line=<FF>){
		chomp($line);
		@tmp = ();
		@tmp = split(/\t/, $line);
		$asId = shift(@tmp);
		$coordinateString = join("\t", @tmp);
		$coordinateToAsIdHref->{$coordinateString}->{"newAsId"} = $asId;
	}
	close FF;
}

# 将新增的AS输出
open A3SS, ">$newAddedA3ssFile";
print A3SS "ASID\tGeneID\tgeneSymbol\tchr\tstrand\tlongExonStart_0base\tlongExonEnd\tshortES\tshortEE\tflankingES\tflankingEE\n";

open A5SS, ">$newAddedA5ssFile";
print A5SS "ASID\tGeneID\tgeneSymbol\tchr\tstrand\tlongExonStart_0base\tlongExonEnd\tshortES\tshortEE\tflankingES\tflankingEE\n";

open RI, ">$newAddedRiFile";
print RI "ASID\tGeneID\tgeneSymbol\tchr\tstrand\triExonStart_0base\triExonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\n";

open MXE, ">$newAddedMxeFile";
print MXE "ASID\tGeneID\tgeneSymbol\tchr\tstrand\t1stExonStart_0base\t1stExonEnd\t2ndExonStart_0base\t2ndExonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\n";

open SE, ">$newAddedSeFile";
print SE "ASID\tGeneID\tgeneSymbol\tchr\tstrand\texonStart_0base\texonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\n";

# 新增ASId编号到最后统一编号，这里升级为1
my ($newAddAsId, $uniformAsId);
open WW, ">$newAsIdToUniformAsIdTsv";
print WW join("\t", "newAddedAsId", "UniformAsId") . "\n";

my @corrdString = keys(%coordinateToAsId);
my ($newAddedA3ssNum, $newAddedA5ssNum, $newAddedRiNum, $newAddedSeNum, $newAddedMxeNum) = (1000000000, 1000000000, 1000000000, 1000000000, 1000000000);
foreach my $coordinateString(@corrdString){
	
	# 如果在新增样本中没有检测到，那么直接放弃
	next if(not(exists($coordinateToAsIdHref->{$coordinateString}->{"newAsId"})));

	# 如果在新样本中检测到的AS
	if($coordinateToAsIdHref->{$coordinateString}->{"newAsId"}=~/A3SS\d+/){

		if(not(exists($coordinateToAsIdHref->{$coordinateString}->{"oldAsId"}))){
			# 新增的AS没有在原来的AS中出现，那么重新编号
			$newAddedA3ssNum++;
			$uniformAsId = $abbr . "A3SS" . $newAddedA3ssNum;
			print A3SS join("\t", $uniformAsId, $coordinateString) . "\n";

		}else{
			# 新增的AS在原来的AS中出现，采用原来的编号
			$uniformAsId = $coordinateToAsIdHref->{$coordinateString}->{"oldAsId"}
		}

		print WW join("\t", $uniformAsId, $coordinateToAsIdHref->{$coordinateString}->{"newAsId"}) . "\n";
		$toUniformAsId{$coordinateToAsIdHref->{$coordinateString}->{"newAsId"}} = $uniformAsId;

	}elsif($coordinateToAsIdHref->{$coordinateString}->{"newAsId"}=~/A5SS\d+/){
			
		if(not(exists($coordinateToAsIdHref->{$coordinateString}->{"oldAsId"}))){
			# 新增的AS没有在原来的AS中出现，那么重新编号
			$newAddedA5ssNum++;
			$uniformAsId = $abbr . "A5SS" . $newAddedA5ssNum;
			print A5SS join("\t", $uniformAsId, $coordinateString) . "\n";

		}else{
			# 新增的AS在原来的AS中出现，采用原来的编号
			$uniformAsId = $coordinateToAsIdHref->{$coordinateString}->{"oldAsId"}
		}

		print WW join("\t", $uniformAsId, $coordinateToAsIdHref->{$coordinateString}->{"newAsId"}) . "\n";
		$toUniformAsId{$coordinateToAsIdHref->{$coordinateString}->{"newAsId"}} = $uniformAsId;

	}elsif($coordinateToAsIdHref->{$coordinateString}->{"newAsId"}=~/MXE\d+/){
			
		if(not(exists($coordinateToAsIdHref->{$coordinateString}->{"oldAsId"}))){
			# 新增的AS没有在原来的AS中出现，那么重新编号
			$newAddedMxeNum++;
			$uniformAsId = $abbr . "MXE" . $newAddedMxeNum;
			print MXE join("\t", $uniformAsId, $coordinateString) . "\n";

		}else{
			# 新增的AS在原来的AS中出现，采用原来的编号
			$uniformAsId = $coordinateToAsIdHref->{$coordinateString}->{"oldAsId"}
		}

		print WW join("\t", $uniformAsId, $coordinateToAsIdHref->{$coordinateString}->{"newAsId"}) . "\n";
		$toUniformAsId{$coordinateToAsIdHref->{$coordinateString}->{"newAsId"}} = $uniformAsId;

	}elsif($coordinateToAsIdHref->{$coordinateString}->{"newAsId"}=~/RI\d+/){
			
		if(not(exists($coordinateToAsIdHref->{$coordinateString}->{"oldAsId"}))){
			# 新增的AS没有在原来的AS中出现，那么重新编号
			$newAddedRiNum++;
			$uniformAsId = $abbr . "RI" . $newAddedRiNum;
			print RI join("\t", $uniformAsId, $coordinateString) . "\n";

		}else{
			# 新增的AS在原来的AS中出现，采用原来的编号
			$uniformAsId = $coordinateToAsIdHref->{$coordinateString}->{"oldAsId"}
		}


		print WW join("\t", $uniformAsId, $coordinateToAsIdHref->{$coordinateString}->{"newAsId"}) . "\n";
		$toUniformAsId{$coordinateToAsIdHref->{$coordinateString}->{"newAsId"}} = $uniformAsId;

	}elsif($coordinateToAsIdHref->{$coordinateString}->{"newAsId"}=~/SE\d+/){
			
		if(not(exists($coordinateToAsIdHref->{$coordinateString}->{"oldAsId"}))){
			# 新增的AS没有在原来的AS中出现，那么重新编号
			$newAddedSeNum++;
			$uniformAsId = $abbr . "SE" . $newAddedSeNum;
			print SE join("\t", $uniformAsId, $coordinateString) . "\n";

		}else{
			# 新增的AS在原来的AS中出现，采用原来的编号
			$uniformAsId = $coordinateToAsIdHref->{$coordinateString}->{"oldAsId"}
		}

		print WW join("\t", $uniformAsId, $coordinateToAsIdHref->{$coordinateString}->{"newAsId"}) . "\n";
		$toUniformAsId{$coordinateToAsIdHref->{$coordinateString}->{"newAsId"}} = $uniformAsId;
	}
}

close WW;
close A3SS;
close A5SS;
close RI;
close MXE;
close SE;

# 打开jcec文件，转换其中的AsId
my (@jcecFile, $jcecFile, $newJcecFile, $line, @tmp);
@jcecFile = split(/,/, $jcecFileList);
foreach $jcecFile(@jcecFile){
	$newJcecFile = $jcecFile . ".with.uniformAs.tsv";
	open WW, ">$newJcecFile";
	open FF, "<$jcecFile";
	$line = <FF>;
	print WW $line;
	while($line=<FF>){
		@tmp = ();
		@tmp = split(/\t/, $line);
		$tmp[0] = $toUniformAsId{$tmp[0]};
		print WW join("\t", @tmp);
	}
	close FF;
	close FF;
}

# 打开jc文件，转换其中的AsId
my (@jcFile, $jcFile, $newJcFile, $line, @tmp);
@jcFile = split(/,/, $jcFileList);
foreach $jcFile(@jcFile){
	$newJcFile = $jcFile . ".with.uniformAs.tsv";
	open WW, ">$newJcFile";
	open FF, "<$jcFile";
	$line = <FF>;
	print WW $line;
	while($line=<FF>){
		@tmp = ();
		@tmp = split(/\t/, $line);
		$tmp[0] = $toUniformAsId{$tmp[0]};
		print WW join("\t", @tmp);
	}
	close FF;
	close FF;
}

