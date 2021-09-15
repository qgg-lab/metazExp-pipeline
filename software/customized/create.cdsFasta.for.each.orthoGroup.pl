#!/usr/bin/perl
use strict;
use Getopt::Long;
use DBI;
if($#ARGV<0){
	print "\nperl $0 \\\n" . 
                "--inputCdsFileList \\\n" .
                "--orthoGroupListFile \\\n" .
		"--pepAlignmentDir \\\n" .
                "--cdsSeqDir \n";
	exit;
}

my ($inputCdsFileList, $orthoGroupListFile, $pepAlignmentDir, $cdsSeqDir);

GetOptions(
        'inputCdsFileList=s'=>\$inputCdsFileList,
        'orthoGroupListFile=s'=>\$orthoGroupListFile,
        'pepAlignmentDir=s'=>\$pepAlignmentDir,
        'cdsSeqDir=s'=>\$cdsSeqDir,
);

# 将所有的CDS序列读入hash
my (@cdsFile, $cdsFile);
my (%cdsSeq, @id, $line, $id);
@cdsFile=split(/,/, $inputCdsFileList);
foreach $cdsFile(@cdsFile){
	open FF, "<$cdsFile";
	while($line=<FF>){
		chomp($line);
		if($line=~/>(.*)/){
			$id = $1;
			@id = split(/ /, $id);
			$id = $id[0];
		}else{
			$cdsSeq{$id}.=$line;
		}
	}
	close FF;
}

print "finish load cds into hash.\n";

# 读取orthoGroup.tsv中的og编号
my ($cmd, $ogIdList, @ogId,  $ogId, $seqIdList, @seqId, $seqId);
$cmd = "grep -v \"Orthogroup\" " . $orthoGroupListFile . " | awk -F \'\\t\' \'{print \$1}\'";
#print $cmd;
#<STDIN>;
$ogIdList = `$cmd`;
@ogId = split(/\n/, $ogIdList);
foreach $ogId(@ogId){
	# 读取pepAlignment文件，获得seqId到数组中
	$cmd = "grep \">\" " . $pepAlignmentDir . "/" . $ogId . ".align | awk -F \'>\' \'{print \$2}\'";
#	print $cmd;
#	<STDIN>;
	$seqIdList = `$cmd`;
	@seqId = ();
	@seqId = split(/\n/, $seqIdList);
	# 按照pepAlingment中序列次序，将对应的CDS序列写入文件
	open WW, ">" . $cdsSeqDir . "/" . $ogId . ".cds.faa";
	foreach $seqId(@seqId){
		print WW ">$seqId\n";
		print WW $cdsSeq{$seqId} . "\n";
	}
	close WW;
}
