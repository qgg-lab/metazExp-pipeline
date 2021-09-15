#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
		"--rmats /mnt/home/liujind1/software/rMATS.4.0.2/rMATS-turbo-Linux-UCS2/rmats.py \\\n" .
                "--subSeqIdListDir ../../006-form-final-trspt-annotation/\\\n" .
                "--subSeqIdListPrefix sub.seqIdList\\\n" .
                "--subGtfDir ../../006-form-final-trspt-annotation/\\\n" .
                "--subGtfPrefix sub.gtf\\\n" .
                "--threadNum 16 \\\n" .
                "--inputBamFile DRX081307.bam\\\n" .
                "--outputDir \\\n" .
                "--outputSubPrefix sub\\\n" .
		"--seqInfoFile\n";
	exit;
}

my (
$cmd,
$samtools, $python, $rmats, $subSeqIdListDir, 
$subSeqIdListPrefix, $subGtfDir, $subGtfPrefix, 
$threadNum, $inputBamFile, $outputDir, $outputSubPrefix,
$seqInfoFile,
$seqLayout, $libraryType, $seqReadLen
);

GetOptions(
	'rmats=s'=>\$rmats,
        'subSeqIdListDir=s'=>\$subSeqIdListDir,
        'subSeqIdListPrefix=s'=>\$subSeqIdListPrefix,
	'subGtfDir=s'=>\$subGtfDir,
	'subGtfPrefix=s'=>\$subGtfPrefix,
	'threadNum=s'=>\$threadNum,
	'inputBamFile=s'=>\$inputBamFile,
	'outputDir=s'=>\$outputDir,
	'outputSubPrefix=s'=>\$outputSubPrefix,
	'seqInfoFile=s'=>\$seqInfoFile,
);

# 排序bam
my $srtBamFile = $inputBamFile . ".srt";
$cmd = "samtools sort -@ " . $threadNum . " -O BAM -o $srtBamFile $inputBamFile 2> $outputDir/log.txt";
system($cmd);


# index bam
$cmd = "samtools index -c -m 14 $srtBamFile 2> $outputDir/log.txt";
system($cmd);


# 获得该experiment的测序信息
&getSeqInfo($seqInfoFile, \$seqReadLen, \$seqLayout, \$libraryType);


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
open FROMGTFMX, ">$fromGtfMxe";
print FROMGTFMX "ID\tGeneID\tgeneSymbol\tchr\tstrand\t1stExonStart_0base\t1stExonEnd\t2ndExonStart_0base\t2ndExonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\n";
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
$subNum = 0;

foreach $subSeqIdListFile(@subSeqIdListFile){

	#print "seqIdList:" . $subSeqIdListFile . "\n";
	#<STDIN>;
	# 获得序号
	my @tt = ();
	@tt = split(/\./, $subSeqIdListFile);
	$subNum = $tt[$#tt];
	# 建立sub输出目录
	$outputSubDir = $outputDir . "/" . $outputSubPrefix . $subNum;
	$cmd = "mkdir -p " . $outputSubDir;
	system($cmd);

	# 输出头
	$cmd = "samtools view -H $srtBamFile > $outputSubDir/out.sam 2> $outputSubDir/log.txt";
	system($cmd);

	# 从seqIdList文件中提取seqId
	$cmd = "cat " . $subSeqIdListFile;
	$seqIdList = `$cmd`;
	@seqId = ();
	@seqId = split(/\n/, $seqIdList);

	# 从srtBam文件中逐个提取比对信息输入
	foreach $seqId(@seqId){
		$cmd = "samtools view $srtBamFile $seqId >> $outputSubDir/out.sam 2> $outputSubDir/log.txt";
		system($cmd);
	}

	# 生成bam文件
	$cmd = "samtools view -b -o $outputSubDir/out.bam $outputSubDir/out.sam 2> $outputSubDir/log.txt";
	system($cmd);

	# 生成b1.txt文件
	my $b1 = "$outputSubDir/b1";
	$cmd = "echo \"$outputSubDir/out.bam\" > $b1";
	system($cmd);
	
	# 执行rmats
	$subGtfFile = $subGtfDir . "/" . $subGtfPrefix . "." . $subNum;
	#print "gtf:" . $subGtfFile . "\n";
	#<STDIN>;
	$cmd = "python $rmats --b1 " . $b1 . " --gtf $subGtfFile --od $outputSubDir -t $seqLayout --nthread $threadNum --tstat $threadNum --libType $libraryType --readLength $seqReadLen --statoff &> $outputSubDir/log.txt";
	system($cmd);

	# 删除bam和sam文件
	$cmd = "rm -rf $outputSubDir/out.bam";
#	system($cmd);
#	$cmd = "rm -rf $outputSubDir/out.sam";
	system($cmd);

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

#	$cmd = "touch " . $outputSubDir . "/rMats.finished";
#	system($cmd);

	#$subNum++;
	#last if($subNum==5);
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

$cmd = "touch " . $outputDir . "/rMats.finished";
system($cmd);


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

