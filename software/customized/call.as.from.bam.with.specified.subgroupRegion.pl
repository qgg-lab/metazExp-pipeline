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
		"--specifiedSubGroupRegion 0-24 \\\n" .
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
$subSeqIdListPrefix, $subGtfDir, $subGtfPrefix, $specifiedSubGroupRegion,
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
	'specifiedSubGroupRegion=s'=>\$specifiedSubGroupRegion,
	'threadNum=s'=>\$threadNum,
	'inputBamFile=s'=>\$inputBamFile,
	'outputDir=s'=>\$outputDir,
	'outputSubPrefix=s'=>\$outputSubPrefix,
	'seqInfoFile=s'=>\$seqInfoFile,
);

# 排序bam
my $srtBamFile = $inputBamFile . ".srt";
if(not -e $inputBamFile . ".srt"){
	$cmd = "samtools sort -@ " . $threadNum . " -O BAM -o $srtBamFile $inputBamFile 2> $outputDir/log.txt";
	system($cmd);
}

# index bam
if(not -e $inputBamFile . ".srt.bai"){
	$cmd = "samtools index -c -m 14 $srtBamFile 2> $outputDir/log.txt";
	system($cmd);
}


# 获得该experiment的测序信息
&getSeqInfo($seqInfoFile, \$seqReadLen, \$seqLayout, \$libraryType);


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

my ($subGroupRegionStart, $subgroupRegionStop)=split(/\-/, $specifiedSubGroupRegion);;
for($subNum=$subGroupRegionStart; $subNum<=$subgroupRegionStop; $subNum++){
#foreach $subSeqIdListFile(@subSeqIdListFile){

	$subSeqIdListFile = $subSeqIdListDir . "/" . $subSeqIdListPrefix . "." . $subNum;

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
	system($cmd);
	$cmd = "rm -rf $outputSubDir/out.sam";
	system($cmd);

	# 标注当前的结束
	$cmd = "touch " . $outputSubDir . "/rMats.finished";
	system($cmd);

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

