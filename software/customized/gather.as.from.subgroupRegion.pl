#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--subSeqIdListDir ../../006-form-final-trspt-annotation/\\\n" .
                "--subSeqIdListPrefix sub.seqIdList\\\n" .
                "--subGtfDir ../../006-form-final-trspt-annotation/\\\n" .
                "--subGtfPrefix sub.gtf\\\n" .
                "--outputDir \\\n" .
                "--outputSubPrefix sub\\\n" .
		"--subGroupRegion 0-49 \\\n" .
		"--seqInfoFile\n";
	exit;
}

my (
$cmd,
$samtools, $python, $rmats, $subSeqIdListDir, 
$subSeqIdListPrefix, $subGtfDir, $subGtfPrefix, $subGroupRegion,
$threadNum, $inputBamFile, $outputDir, $outputSubPrefix,
$seqInfoFile,
$seqLayout, $libraryType, $seqReadLen
);

GetOptions(
        'subSeqIdListDir=s'=>\$subSeqIdListDir,
        'subSeqIdListPrefix=s'=>\$subSeqIdListPrefix,
	'subGtfDir=s'=>\$subGtfDir,
	'subGtfPrefix=s'=>\$subGtfPrefix,
	'outputDir=s'=>\$outputDir,
	'outputSubPrefix=s'=>\$outputSubPrefix,
	'subGroupRegion=s'=>\$subGroupRegion,
	'seqInfoFile=s'=>\$seqInfoFile,
);

# 定义最后输出的4类结果文件
# A5SS
my ($fromGtfA5ss, $fromGtfNovelA5ss, $jcecA5ss, $jcA5ss) = ($outputDir . "/fromGTF.A5SS.txt", $outputDir . "/fromGTF.novelEvents.A5SS.txt", $outputDir . "/JCEC.raw.input.A5SS.txt", $outputDir . "/JC.raw.input.A5SS.txt");
open FROMGTFA5SS, ">$fromGtfA5ss";
print FROMGTFA5SS "ID\tGeneID\tgeneSymbol\tchr\tstrand\tlongExonStart_0base\tlongExonEnd\tshortES\tshortEE\tflankingES\tflankingEE\n";
open FROMGTFNOVELA5SS, ">$fromGtfNovelA5ss";
print FROMGTFNOVELA5SS "ID\tGeneID\tgeneSymbol\tchr\tstrand\tlongExonStart_0base\tlongExonEnd\tshortES\tshortEE\tflankingES\tflankingEE\n";
open JCECA5SS, ">$jcecA5ss";
print JCECA5SS "ID\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\n";
open JCA5SS, ">$jcA5ss";
print JCA5SS "ID\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\n";
# A3SS
my ($fromGtfA3ss, $fromGtfNovelA3ss, $jcecA3ss, $jcA3ss) = ($outputDir . "/fromGTF.A3SS.txt", $outputDir . "/fromGTF.novelEvents.A3SS.txt", $outputDir . "/JCEC.raw.input.A3SS.txt", $outputDir . "/JC.raw.input.A3SS.txt");
open FROMGTFA3SS, ">$fromGtfA3ss";
print FROMGTFA3SS "ID\tGeneID\tgeneSymbol\tchr\tstrand\tlongExonStart_0base\tlongExonEnd\tshortES\tshortEE\tflankingES\tflankingEE\n";
open FROMGTFNOVELA3SS, ">$fromGtfNovelA3ss";
print FROMGTFNOVELA3SS "ID\tGeneID\tgeneSymbol\tchr\tstrand\tlongExonStart_0base\tlongExonEnd\tshortES\tshortEE\tflankingES\tflankingEE\n";
open JCECA3SS, ">$jcecA3ss";
print JCECA3SS "ID\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\n";
open JCA3SS, ">$jcA3ss";
print JCA3SS "ID\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\n";
# MXE
my ($fromGtfMxe, $fromGtfNovelMxe, $jcecMxe, $jcMxe) = ($outputDir . "/fromGTF.MXE.txt", $outputDir . "/fromGTF.novelEvents.MXE.txt", $outputDir . "/JCEC.raw.input.MXE.txt", $outputDir . "/JC.raw.input.MXE.txt");
open FROMGTFMXE, ">$fromGtfMxe";
print FROMGTFMXE "ID\tGeneID\tgeneSymbol\tchr\tstrand\t1stExonStart_0base\t1stExonEnd\t2ndExonStart_0base\t2ndExonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\n";
open FROMGTFNOVELMXE, ">$fromGtfNovelMxe";
print FROMGTFNOVELMXE "ID\tGeneID\tgeneSymbol\tchr\tstrand\t1stExonStart_0base\t1stExonEnd\t2ndExonStart_0base\t2ndExonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\n";
open JCECMXE, ">$jcecMxe";
print JCECMXE "ID\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\n";
open JCMXE, ">$jcMxe";
print JCMXE "ID\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\n";

# RI
my ($fromGtfRi, $fromGtfNovelRi, $jcecRi, $jcRi) = ($outputDir . "/fromGTF.RI.txt", $outputDir . "/fromGTF.novelEvents.RI.txt", $outputDir . "/JCEC.raw.input.RI.txt", $outputDir . "/JC.raw.input.RI.txt");
open FROMGTFRI, ">$fromGtfRi";
print FROMGTFRI "ID\tGeneID\tgeneSymbol\tchr\tstrand\triExonStart_0base\triExonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\n";
open FROMGTFNOVELRI, ">$fromGtfNovelRi";
print FROMGTFNOVELRI "ID\tGeneID\tgeneSymbol\tchr\tstrand\triExonStart_0base\triExonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\n";
open JCECRI, ">$jcecRi";
print JCECRI "ID\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\n";
open JCRI, ">$jcRi";
print JCRI "ID\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\n";

# SE
my ($fromGtfSe, $fromGtfNovelSe, $jcecSe, $jcSe) = ($outputDir . "/fromGTF.SE.txt", $outputDir . "/fromGTF.novelEvents.SE.txt", $outputDir . "/JCEC.raw.input.SE.txt", $outputDir . "/JC.raw.input.SE.txt");
open FROMGTFSE, ">$fromGtfSe";
print FROMGTFSE "ID\tGeneID\tgeneSymbol\tchr\tstrand\texonStart_0base\texonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\n";
open FROMGTFNOVELSE, ">$fromGtfNovelSe";
print FROMGTFNOVELSE "ID\tGeneID\tgeneSymbol\tchr\tstrand\texonStart_0base\texonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\n";
open JCECSE, ">$jcecSe";
print JCECSE "ID\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\n";
open JCSE, ">$jcSe";
print JCSE "ID\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\n";



###########################################################################
# 获得当前基因组分组后的所有subSeqIdList文件
my ($subNum, $subSeqIdListFileList, @subSeqIdListFile, $subSeqIdListFile);
my ($subGtfFile);
my ($seqIdList, @seqId, $seqId, $line, $asId, @field);
my ($outputSubDir);
# ../../006-form-final-trspt-annotation/sub.seqIdList.0
# ../../006-form-final-trspt-annotation/sub.seqIdList.1
# ../../006-form-final-trspt-annotation/sub.seqIdList.10
# ../../006-form-final-trspt-annotation/sub.seqIdList.11
$cmd = "ls -1 " . $subSeqIdListDir . "/" . $subSeqIdListPrefix . "*";
$subSeqIdListFileList = `$cmd`;
@subSeqIdListFile=split(/\n/, $subSeqIdListFileList);

my ($subGroupRegionStart, $subGroupRegionStop) = split(/\-/, $subGroupRegion);

$subNum = 0;
my $allSuccessFlag = 1;
for($subNum=$subGroupRegionStart; $subNum<=$subGroupRegionStop; $subNum++){

	# 建立sub目录是否存在
	$outputSubDir = $outputDir . "/" . $outputSubPrefix . $subNum;
	if(not -e $outputSubDir){
		$allSuccessFlag = 0;
		last;
	}

	# 检测不到rMats.finished，那么就提前结束
	if(not -e $outputSubDir . "/rMats.finished"){
		$allSuccessFlag = 0;
		last;
	}

	# 将4类文件：
	# fromGTF.A5SS, fromGTF.novelEvents.A5SS.txt, JCEC.raw.input.A5SS.txt, JC.raw.input.A5SS.txt
	# 提取出来输出到分别输出到一个文件, 注意将每行的编号改变一下，规则为：原来Id . "_". subNum
	# fromGTF.A5SS.txt
	open FF, "<$outputSubDir/fromGTF.A5SS.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print FROMGTFA5SS join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;
	open FF, "<$outputSubDir/fromGTF.novelEvents.A5SS.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print FROMGTFNOVELA5SS join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;
	open FF, "<$outputSubDir/JCEC.raw.input.A5SS.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print JCECA5SS join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;
	open FF, "<$outputSubDir/JC.raw.input.A5SS.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print JCA5SS join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;
	# fromGTF.A3SS.txt
	open FF, "<$outputSubDir/fromGTF.A3SS.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print FROMGTFA3SS join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;
	open FF, "<$outputSubDir/fromGTF.novelEvents.A3SS.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print FROMGTFNOVELA3SS join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;
	open FF, "<$outputSubDir/JCEC.raw.input.A3SS.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print JCECA3SS join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;
	open FF, "<$outputSubDir/JC.raw.input.A3SS.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print JCA3SS join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;
	# fromGTF.MXE.txt
	open FF, "<$outputSubDir/fromGTF.MXE.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print FROMGTFMXE join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;
	open FF, "<$outputSubDir/fromGTF.novelEvents.MXE.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print FROMGTFNOVELMXE join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;
	open FF, "<$outputSubDir/JCEC.raw.input.MXE.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print JCECMXE join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;
	open FF, "<$outputSubDir/JC.raw.input.MXE.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print JCMXE join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;
	# fromGTF.RI.txt
	open FF, "<$outputSubDir/fromGTF.RI.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print FROMGTFRI join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;
	open FF, "<$outputSubDir/fromGTF.novelEvents.RI.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print FROMGTFNOVELRI join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;
	open FF, "<$outputSubDir/JCEC.raw.input.RI.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print JCECRI join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;
	open FF, "<$outputSubDir/JC.raw.input.RI.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print JCRI join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;
	# fromGTF.SE.txt
	open FF, "<$outputSubDir/fromGTF.SE.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print FROMGTFSE join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;
	open FF, "<$outputSubDir/fromGTF.novelEvents.SE.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print FROMGTFNOVELSE join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;
	open FF, "<$outputSubDir/JCEC.raw.input.SE.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print JCECSE join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;
	open FF, "<$outputSubDir/JC.raw.input.SE.txt";
	<FF>;
	while($line=<FF>){
		@field = ();
		@field = split(/\t/, $line);
		$asId = shift(@field);
		print JCSE join("\t", $asId . "_" . $subNum, @field);
	}
	close FF;

	$cmd = "rm -rf " . $outputSubDir;
#	system($cmd);

}

# 关闭文件
close FROMGTFA5SS;
close FROMGTFNOVELA5SS;
close JCECA5SS;
close JCA5SS;
close FROMGTFA3SS;
close FROMGTFNOVELA3SS;
close JCECA3SS;
close JCA3SS;
close FROMGTFMXE;
close FROMGTFNOVELMXE;
close JCECMXE;
close JCMXE;
close FROMGTFRI;
close FROMGTFNOVELRI;
close JCECRI;
close JCRI;
close FROMGTFSE;
close FROMGTFNOVELSE;
close JCECSE;
close JCSE;

# 检测每个sub是否全部正常结束
if($allSuccessFlag == 1){
	$cmd = "touch " . $outputDir . "/rMats.finished";
	system($cmd);
}else{
	system("rm -rf " . $outputDir . "/rMats.finished");
	system("rm -rf " . $outputDir . "/fromGTF.*");
	system("rm -rf " . $outputDir . "/JCEC.raw.*");
	system("rm -rf " . $outputDir . "/JC.raw.*");
}


sub getSeqInfo{
        my ($alignmentFile, $seqReadLen, $seqLayout, $libraryType) = @_;
        my ($line, @field);
        open AA, "<$alignmentFile";
        <AA>;
        $line=<AA>;
        chomp($line);
        @field = split(/\t/, $line);
        $$seqReadLen = $field[3];
        $$seqLayout = lc($field[2]);
        $$libraryType = $field[6];
        if($$libraryType eq "RF" or $$libraryType eq "R"){
                $$libraryType = "fr-firststrand";
        }elsif($$libraryType eq "FR" or $$libraryType eq "F"){
                $$libraryType = "fr-secondstrand";
        }else{
                $$libraryType = "fr-unstranded";
        }
}

