#!/usr/bin/perl
use strict;
use Getopt::Long;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--inputAsFileList whole.A5SS.catalog,whole.SE.catalog,whole.MXE.catalog,whole.RI.catalog,whole.A3SS.catalog\\\n" .
                "--inputJcecFileList whole.jcec.A5SS.tsv,whole.jcec.SE.tsv,whole.jcec.MXE.tsv,whole.jcec.RI.tsv,whole.jcec.A3SS.tsv\\\n" .
                "--inputSampleListFile selected.experiment.list.tsv\\\n" .
                "--minIncluReadNumInOneExp 3 \\\n" .
                "--minSkipingReadNumInOneExp 3 \\\n" .
                "--minPrjNum 1 \\\n" .
                "--minExpNum 3 \\\n" .
		"--outputJcecFileList jcec.A5SS.tsv,jcec.SE.tsv,jcec.MXE.tsv,jcec.RI.tsv,jcec.A3SS.tsv \\\n" .
                "--outputAsFileList A5SS.catalog,SE.catalog,MXE.catalog,RI.catalog,A3SS.catalog\\\n" .
		"--workDir ./\n";
	exit;
}

my ($inputAsFileList, $inputJcecFileList, $inputSampleListFile, $minIncluReadNumInOneExp);
my ($minSkipingReadNumInOneExp, $minPrjNum, $minExpNum, $outputAsFileList, $outputJcecFileList, $workDir);
GetOptions(
	'inputAsFileList=s'=>\$inputAsFileList, 
	'inputJcecFileList=s'=>\$inputJcecFileList,
	'inputSampleListFile=s'=>\$inputSampleListFile,
	'minIncluReadNumInOneExp=s'=>\$minIncluReadNumInOneExp,
	'minSkipingReadNumInOneExp=s'=>\$minSkipingReadNumInOneExp,
	'minPrjNum=s'=>\$minPrjNum,
	'minExpNum=s'=>\$minExpNum,
	'outputAsFileList=s'=>\$outputAsFileList,
	'outputJcecFileList=s'=>\$outputJcecFileList,
	'workDir=s'=>\$workDir,
);

system("mkdir -p " . $workDir);
my (@wholeAsFile, $wholeAsFile, @jcecFile, $jcecFile, @outputAsFile, $outputAsFile, @outputJcecFile, $outputJcecFile);
my (@field, $field, $line, $i);
my (%expToPrj);
open FF, "<$inputSampleListFile";
while($line=<FF>){
	#9031    SRP097223       SRX2506291      33      RF      PAIRED  125     24.21   95.69%  23.166549
	@field = split(/\t/, $line);
	$expToPrj{$field[2]}=$field[1];
}
close FF;

my (%as, @asId, $asId, $inclusionNum, $skippingNum);
my ($expId, $psi);
@wholeAsFile = split(/,/, $inputAsFileList);
@jcecFile = split(/,/, $inputJcecFileList);
@outputAsFile = split(/,/, $outputAsFileList);
@outputJcecFile = split(/,/, $outputJcecFileList);

for($i=0; $i<=$#wholeAsFile; $i++){
	$wholeAsFile = $wholeAsFile[$i];
	$jcecFile = $jcecFile[$i];
	$outputAsFile = $outputAsFile[$i];
	$outputJcecFile = $outputJcecFile[$i];

	%as = ();
	# read jcec
	open FF, "<$jcecFile";
	open JCEC, ">$outputJcecFile";
	$line=<FF>;
	chomp($line);
	print JCEC $line . "\tPsi\n";
	#ASID    IJC_SAMPLE_1    SJC_SAMPLE_1    IncFormLen      SkipFormLen     Experiment  PSI
	while($line=<FF>){
		#GGALA5SS0000000001      527     31      730     100     SRX1036607
		chomp($line);
		@field = ();
		@field = split(/\t/, $line);
		$asId = $field[0];
		$inclusionNum = $field[1];
		$skippingNum = $field[2];
		$expId = $field[5];
		if(($inclusionNum >= $minIncluReadNumInOneExp or $skippingNum >= $minSkipingReadNumInOneExp) and exists($expToPrj{$expId})){
			if($field[1]==0 and $field[2] !=0){
				$psi = 0;
			}elsif($field[1]==0 and $field[2] ==0){
				$psi = "NA";
			}elsif($field[1]!=0 and $field[2] ==0){
				$psi = 1;
			}else{
				$psi = $field[1]/$field[3]/($field[1]/$field[3]+$field[2]/$field[4]);
			}
			
			print JCEC $line . "\t" . $psi . "\n";
			${$as{$asId}}{"expNum"}++;
			if(not(exists(${$as{$asId}}{"prjList"}))){
				${$as{$asId}}{"prjList"} = $expToPrj{$expId};
				${$as{$asId}}{"prjNum"} = 1;
			}else{
				if(index(${$as{$asId}}{"prjList"}, $expToPrj{$expId})<0){
					${$as{$asId}}{"prjList"} .= "#" . $expToPrj{$expId};
					${$as{$asId}}{"prjNum"}++;
				}
			}
		}
	}
	close FF;
	close JCEC;

	##########################################
	# 读取AS，如果存在于as哈希中，那么输出
	open FF, "<$wholeAsFile";
	open WW, ">$outputAsFile";
	open WID, ">$workDir/asId.txt";
	#ASID    GeneID  geneSymbol      chr     strand  longExonStart_0base     longExonEnd     shortES shortEE flankingES      flankingEE
	$line = <FF>;
	print WW $line;
	while($line=<FF>){
		#GGALA5SS0000014947 "ENSGALG00000000003" "PANX2" chr1 + 20329069 20345783 20329069 20329191 20347539 20349394
		@field = split(/\t/, $line);
		$asId = $field[0];
		if(${$as{$asId}}{"expNum"} >= $minExpNum and ${$as{$asId}}{"prjNum"} >= $minPrjNum){
			print WW $line;
			print WID $asId . "\n";
		}
	}
	close WW;
	close FF;
	close WID;

	# 在jcec文件删除不在输出AS文件中jcec
	system("grep ASID " . $outputJcecFile . " > " . $workDir . "/tmp.jcec.tsv");
	system("grep -Ff " . $workDir . "/asId.txt " . $outputJcecFile . " >>" . $workDir . "/tmp.jcec.tsv");
	system("rm -rf $outputJcecFile");
	system("mv " . $workDir . "/tmp.jcec.tsv " . $outputJcecFile);
}
